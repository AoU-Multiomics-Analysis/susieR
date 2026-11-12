// Use the Cromwell 85 command renderer from the pinned WOMtool JAR.
// Do not start a workflow, Docker, PLINK, model fits, or a Terra run.
import java.nio.charset.StandardCharsets.UTF_8
import java.nio.file.{Files, Paths}
import com.typesafe.config.ConfigFactory
import languages.wdl.draft3.WdlDraft3LanguageFactory
import scala.collection.mutable
import scala.jdk.CollectionConverters._
import scala.sys.process._
import scala.concurrent.Future
import wom.callable.{CommandTaskDefinition, RuntimeEnvironment}
import eu.timepit.refined.types.numeric.PosInt
import wom.expression.EmptyIoFunctionSet
import wom.types._
import wom.values._
import wom.WomFileMapper

object FullAggregateWomRender extends App {
  val root = Paths.get(args.headOption.getOrElse(".")).toAbsolutePath.normalize()
  val factory = new WdlDraft3LanguageFactory(ConfigFactory.parseString("enabled = true"))
  val descriptors = Vector(
    "workflows/SusieR2.wdl" -> "AggregatePgenSusie",
    "workflows/AggregateFullSusie.wdl" -> "AggregateFullSusieFiles")
  for ((descriptor, taskName) <- descriptors) {
    val source = Files.readString(root.resolve(descriptor))
    val bundle = factory.getWomBundle(source, None, "", List.empty, List(factory))
      .fold(errors => sys.error(errors.toString), identity)
    val task = bundle.allCallables(taskName).asInstanceOf[CommandTaskDefinition]
    val cwd = Files.createTempDirectory("susier-full-aggregate-wom-")
    try {
      val mapping = mutable.Map.empty[String, String]
      val expectedFiles = mutable.Map.empty[String, String]
      val expectedLists = mutable.Map.empty[String, Vector[String]]
      val expectedMetadata = mutable.Map.empty[String, Vector[String]]
      val generatedFiles = mutable.Set.empty[String]
      val mappedFiles = mutable.Set.empty[String]
      def inputValue(typ: WomType, name: String): WomValue = typ match {
        case WomSingleFileType =>
          val uri = "gs://bucket/inputs/" + name
          val local = cwd.resolve("input 'quoted' \"double\" $(touch PWNED) `touch PWNED` " + name)
          Files.writeString(local, "fixture\n")
          mapping(uri) = local.toString
          expectedFiles(name) = local.toString
          WomSingleFile(uri)
        case array: WomArrayType if array.memberType == WomSingleFileType =>
          val members = Vector(0, 1).map(index => inputValue(array.memberType, name + "_" + index))
          expectedLists(name) = Vector(0, 1).map(index => expectedFiles.remove(name + "_" + index).get)
          WomArray(array, members)
        case array: WomArrayType if array.memberType == WomStringType =>
          val values = if (name == "pgen_chromosomes") Vector("1", "X") else Vector(
            "gs://bucket/outputs/" + name + "/first",
            "gs://bucket/outputs/" + name + "/second 'quoted' \"double\" $(touch PWNED) `touch PWNED`")
          expectedMetadata(name) = values
          WomArray(array, values.map(WomString(_)))
        case WomStringType => WomString("image 'quoted' $(touch PWNED) `touch PWNED`")
        case WomIntegerType => WomInteger(if (name == "memory_gb") 8 else if (name == "disk_gb") 50 else 1)
        case other => sys.error("Unsupported aggregate fixture input type: " + other)
      }
      val inputs = task.inputs.map(input => input -> inputValue(input.womType, input.name)).toMap
      val fullOnly = taskName == "AggregateFullSusieFiles"
      val fileLists = if (fullOnly) Set("full_susie_files") else
        Set("summary_files", "credible_set_files", "index_files", "diagnostic_files")
      assert(expectedLists.keySet.toSet == fileLists,
        "Unexpected localized File arrays for " + taskName + ": " + expectedLists.keySet)
      val metadataUris = expectedMetadata.valuesIterator.flatten.filter(_.startsWith("gs://")).toSet
      // Keep gs:// values in the input environment. The real renderer must
      // map File values before sub() and command-time write_lines run.
      def localize(value: WomValue): WomValue = WomFileMapper.mapWomFiles(file => {
        assert(!metadataUris.contains(file.value), "Manifest URI metadata was treated as a File: " + file.value)
        mappedFiles += file.value
        val path = mapping.getOrElse(file.value,
          if (mapping.valuesIterator.contains(file.value) || generatedFiles.contains(file.value)) file.value
          else sys.error("Unexpected command File: " + file.value))
        WomSingleFile(path)
      })(value).get
      val runtime = RuntimeEnvironment(cwd.toString, cwd.toString, PosInt.unsafeFrom(1), 8.0, 0L, 0L)
      val io = new EmptyIoFunctionSet {
        override def writeFile(path: String, content: String): Future[WomSingleFile] = {
          val output = cwd.resolve(path)
          Files.writeString(output, content)
          generatedFiles += output.toString
          Future.successful(WomSingleFile(output.toString))
        }
      }
      val command = task.instantiateCommand(inputs, io, localize, runtime)
        .toEither.fold(errors => sys.error(errors.toString), _.commandString)
      Files.writeString(cwd.resolve("command.sh"), command)
      val bin = Files.createDirectory(cwd.resolve("bin"))
      val stub = bin.resolve("Rscript")
      Files.writeString(stub, """#!/usr/bin/env python3
import os, pathlib, sys
args = sys.argv[1:]
pathlib.Path('capture.args').write_bytes(b''.join(os.fsencode(arg) + b'\0' for arg in args))
list_flags = ('--summary-files', '--full-susie-files', '--credible-set-files', '--index-files', '--diagnostic-files',
              '--gene-rds-uris-file', '--gene-parquet-uris-file', '--gene-tsv-uris-file',
              '--gene-log-uris-file', '--gene-provenance-uris-file', '--pgen-chromosomes-file',
              '--pgen-uris-file', '--pvar-uris-file', '--psam-uris-file')
for flag in list_flags:
    if flag in args:
        source = pathlib.Path(args[args.index(flag) + 1])
        pathlib.Path('capture.' + flag[2:]).write_bytes(source.read_bytes())
""")
      assert(stub.toFile.setExecutable(true), "Could not make the recording CLI executable")
      val stdout = new StringBuilder
      val stderr = new StringBuilder
      val exit = Process(Seq("bash", "command.sh"), cwd.toFile,
        "PATH" -> (bin.toString + java.io.File.pathSeparator + sys.env("PATH")))
        .!(ProcessLogger(line => stdout.append(line).append('\n'), line => stderr.append(line).append('\n')))
      assert(exit == 0, "Cromwell-rendered " + taskName + " command failed:\n" + stderr)
      val actual = new String(Files.readAllBytes(cwd.resolve("capture.args")), UTF_8)
        .split("\u0000", -1).dropRight(1).toVector
      assert(actual.head == "/opt/r/scripts/aggregate_pgen_susie.R", "Unexpected aggregate CLI: " + actual)
      assert(actual.contains("--full-only") == fullOnly, "Incorrect full-only mode")
      assert(actual.contains("--full-susie-files") == fullOnly,
        "Full-result inputs must be used only by the standalone aggregate")
      for ((name, path) <- expectedFiles) {
        val flag = "--" + name.replace('_', '-')
        val index = actual.indexOf(flag)
        assert(index >= 0 && index + 1 < actual.size && actual(index + 1) == path,
          "Cromwell did not pass the exact local path for " + flag + ": " + actual)
        assert(Files.isReadable(Paths.get(path)), "Single File input is not readable: " + path)
      }
      for ((name, paths) <- expectedLists) {
        val flag = "--" + name.replace('_', '-')
        assert(actual.contains(flag), "Missing localized input list: " + flag)
        val captured = Files.readAllLines(cwd.resolve("capture." + flag.drop(2))).asScala.toVector
        assert(captured == paths, "Cromwell did not localize File array members before write_lines: " + captured)
        assert(captured.forall(path => Files.isReadable(Paths.get(path))), "Input list contains an unreadable file")
        assert(!captured.exists(_.contains("://")), "Input list contains an unresolved cloud URI")
      }
      for ((name, uris) <- expectedMetadata) {
        val flag = "--" + name.replace('_', '-') + "-file"
        assert(actual.contains(flag), "Missing output URI metadata list: " + flag)
        val captured = Files.readAllLines(cwd.resolve("capture." + flag.drop(2))).asScala.toVector
        assert(captured == uris, "Cromwell changed output URI metadata: " + captured)
        assert(captured.filter(_.startsWith("gs://")).forall(uri => !mappedFiles.contains(uri)),
          "Output URI metadata was passed through File localization")
      }
      assert(mapping.keySet.subsetOf(mappedFiles), "Some input Files bypassed command localization")
      assert(!actual.exists(_.contains("gs://")), "Unresolved cloud URI reached a CLI argument")
      assert(!Files.exists(cwd.resolve("PWNED")), "File or metadata text executed shell code")
      assert(stdout.toString.contains("COMPLETE"), "The aggregate command did not complete")
      println("Cromwell 85 " + taskName + " command localization and URI metadata checks passed")
    } finally {
      val paths = Files.walk(cwd)
      try paths.iterator().asScala.toVector.sortBy(_.getNameCount).reverse.foreach(Files.delete)
      finally paths.close()
    }
  }
}
