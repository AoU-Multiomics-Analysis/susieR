// Run the real Cromwell 85 command renderer from the pinned WOMtool JAR.
// This checks URI-to-local-path conversion and shell argument quoting. It does
// not start Cromwell, Docker, PLINK, model fits, or a Terra workflow.
import java.nio.charset.StandardCharsets.UTF_8
import java.nio.file.{Files, Paths}
import com.typesafe.config.ConfigFactory
import languages.wdl.draft3.WdlDraft3LanguageFactory
import scala.collection.mutable
import scala.jdk.CollectionConverters._
import scala.sys.process._
import wom.callable.{CommandTaskDefinition, RuntimeEnvironment}
import eu.timepit.refined.types.numeric.PosInt
import wom.expression.EmptyIoFunctionSet
import wom.types._
import wom.values._
import wom.WomFileMapper

object PgenWomRender extends App {
  val root = Paths.get(args.headOption.getOrElse(".")).toAbsolutePath.normalize()
  val source = Files.readString(root.resolve("workflows/SusieR2.wdl"))
  val factory = new WdlDraft3LanguageFactory(ConfigFactory.parseString("enabled = true"))
  val bundle = factory.getWomBundle(source, None, "", List.empty, List(factory))
    .fold(errors => sys.error(errors.toString), identity)
  val task = bundle.allCallables("RunPgenSusieBatch").asInstanceOf[CommandTaskDefinition]

  for (includeOptional <- List(false, true)) {
    val cwd = Files.createTempDirectory("susier-wom-render-")
    try {
      val mapping = mutable.Map.empty[String, String]
      val expectedFiles = mutable.Map.empty[String, String]
      def inputValue(typ: WomType, name: String): WomValue = typ match {
        case WomSingleFileType =>
          val uri = "gs://bucket/" + name
          val local = cwd.resolve("input 'quoted' \"double\" $(touch PWNED) `touch PWNED` " + name)
          Files.writeString(local, "fixture\n")
          mapping(uri) = local.toString
          expectedFiles(name) = local.toString
          WomSingleFile(uri)
        case optional: WomOptionalType =>
          WomOptionalValue(optional.memberType,
            if (includeOptional) Some(inputValue(optional.memberType, name)) else None)
        // Keep support for the original struct input so this test reproduces
        // its failure when run against the old descriptor.
        case composite: WomCompositeType =>
          WomObject(composite.typeMap.map { case (field, fieldType) =>
            field -> inputValue(fieldType, field)
          }, composite)
        case WomBooleanType => WomBoolean(includeOptional)
        case WomIntegerType => WomInteger(if (name == "mac") {
          if (includeOptional) 5 else 0
        } else 1)
        case WomFloatType => WomFloat(0.0)
        case WomStringType =>
          WomString(if (name == "missing_code") "-1" else "qval 'quoted' $(touch PWNED) `touch PWNED`")
        case other => sys.error("Unsupported fixture input type: " + other)
      }
      val inputs = task.inputs.map(input => input -> inputValue(input.womType, input.name)).toMap
      // Do not rewrite the input environment before expression evaluation.
      // Cromwell must map File values while it instantiates the command.
      def localize(value: WomValue): WomValue = WomFileMapper.mapWomFiles(file => {
        val path = mapping.getOrElse(file.value,
          if (mapping.valuesIterator.contains(file.value)) file.value
          else sys.error("Unexpected command File: " + file.value))
        WomSingleFile(path)
      })(value).get
      val runtime = RuntimeEnvironment(cwd.toString, cwd.toString, PosInt.unsafeFrom(1), 16.0, 0L, 0L)
      val command = task.instantiateCommand(inputs, new EmptyIoFunctionSet, localize, runtime)
        .toEither.fold(errors => sys.error(errors.toString), _.commandString)
      Files.writeString(cwd.resolve("command.sh"), command)
      val bin = Files.createDirectory(cwd.resolve("bin"))
      val stub = bin.resolve("Rscript")
      Files.writeString(stub, """#!/usr/bin/env python3
import os, pathlib, sys
pathlib.Path('capture.args').write_bytes(b''.join(os.fsencode(arg) + b'\0' for arg in sys.argv[1:]))
""")
      assert(stub.toFile.setExecutable(true), "Could not make the recording CLI executable")
      val stdout = new StringBuilder
      val stderr = new StringBuilder
      val exit = Process(Seq("bash", "command.sh"), cwd.toFile,
        "PATH" -> (bin.toString + java.io.File.pathSeparator + sys.env("PATH")))
        .!(ProcessLogger(line => stdout.append(line).append('\n'), line => stderr.append(line).append('\n')))
      assert(exit == 0, "Cromwell-rendered batch command failed (optional inputs=" + includeOptional + "):\n" + stderr)
      val actual = new String(Files.readAllBytes(cwd.resolve("capture.args")), UTF_8)
        .split("\u0000", -1).dropRight(1).toVector
      assert(actual.head == "/opt/r/scripts/run_pgen_susie_batch.R", "Unexpected batch CLI")
      for ((name, path) <- expectedFiles) {
        val flag = "--" + name.replace('_', '-')
        val index = actual.indexOf(flag)
        assert(index >= 0 && index + 1 < actual.size && actual(index + 1) == path,
          "Cromwell did not pass the exact local path for " + flag + ": " + actual)
      }
      for (name <- List("sample_list", "ancestry", "variant_list", "additional_genotypes")) {
        assert(actual.contains("--" + name.replace('_', '-')) == includeOptional,
          "Incorrect optional file argument: " + name)
      }
      assert(actual.contains("--missing-code=-1"), "Negative missing code changed")
      assert(actual.contains(if (includeOptional) "--mac=5" else "--mac=0"),
        "MAC threshold did not reach the batch CLI")
      assert(!actual.exists(_.contains("gs://")), "Unresolved cloud input reached the batch CLI")
      assert(!Files.exists(cwd.resolve("PWNED")), "File or string argument executed shell code")
      println("Cromwell 85 batch command localization passed (optional inputs=" + includeOptional + ")")
    } finally {
      val paths = Files.walk(cwd)
      try paths.iterator().asScala.toVector.sortBy(_.getNameCount).reverse.foreach(Files.delete)
      finally paths.close()
    }
  }
}
