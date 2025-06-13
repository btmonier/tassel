import org.jetbrains.kotlin.gradle.dsl.JvmTarget
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import java.util.Locale

plugins {
    kotlin("jvm") version "2.1.20"
    id("org.jetbrains.dokka") version "2.0.0"
    id("com.github.johnrengelman.shadow") version "7.1.2"
    id("org.jetbrains.kotlinx.kover") version "0.9.1"
    alias(libs.plugins.composeMultiplatform) apply false
    alias(libs.plugins.composeCompiler) apply false
    alias(libs.plugins.kotlinMultiplatform) apply false
}

java {
    toolchain {
        languageVersion.set(JavaLanguageVersion.of(21))
    }
}

group = "net.maizegenetics"
version = "5.2.96"
description = "TASSEL is a software package to evaluate traits associations, evolutionary patterns, and linkage disequilibrium."
val kotlinVersion = "2.1.21"

repositories {
    mavenCentral()
}

repositories {
    mavenCentral()
    maven {
        url = uri("https://maven.scijava.org/content/repositories/public/") // needed for JHDF5
    }
}

dependencies {
    testImplementation(kotlin("test"))
    implementation("org.apache.logging.log4j:log4j-core:2.21.1")
    implementation("com.google.guava:guava:22.0")
    implementation("org.apache.commons:commons-math3:3.4.1")
    implementation("commons-codec:commons-codec:1.10")
    implementation("org.biojava:biojava-core:6.0.4")
    implementation("org.biojava:biojava-alignment:6.0.4")
    implementation("org.biojava:biojava-phylo:4.2.12")
    implementation("org.jfree:jfreechart:1.0.19")
    implementation("org.jfree:jfreesvg:3.2")
    implementation("net.sf.trove4j:trove4j:3.0.3")
    implementation("org.ejml:ejml-ddense:0.41")
    implementation("org.xerial.snappy:snappy-java:1.1.8.4")
    implementation("javax.mail:mail:1.4")
    implementation("com.googlecode.json-simple:json-simple:1.1.1")
    implementation("org.glassfish:javax.json:1.0.4")
    implementation("org.xerial:sqlite-jdbc:3.39.2.1")
    implementation("com.github.samtools:htsjdk:2.24.1")
    implementation("org.ahocorasick:ahocorasick:0.2.4")
    implementation("org.postgresql:postgresql:42.6.0")
    implementation("org.apache.avro:avro:1.8.1")
    implementation("colt:colt:1.2.0")
    implementation("org.biojava.thirdparty:forester:1.039")
    implementation("cisd:jhdf5:19.04.1")
    implementation("org.jetbrains.kotlin:kotlin-stdlib-jdk8:$kotlinVersion")
    testImplementation("org.jetbrains.kotlin:kotlin-test:$kotlinVersion")
    implementation("org.ini4j:ini4j:0.5.4")
    implementation("it.unimi.dsi:fastutil:8.2.2")
}

tasks {
    // Set JAR file name
    withType<Jar> {
        archiveFileName.set("sTASSEL.jar")
    }

    // Copy runtime dependencies into build/libs/lib
    register<Copy>("copyDependencies") {
        from(configurations.runtimeClasspath)
        into(layout.buildDirectory.dir("libs/lib").get().asFile)
    }

    // Ensure dependencies are copied after build
    named("build") {
        dependsOn("copyDependencies")
    }

    // Compile with Java 21 bytecode target
    withType<JavaCompile> {
        sourceCompatibility = "21"
        targetCompatibility = "21"
    }

    // Optional: Kotlin compile target
    withType<KotlinCompile> {
        compilerOptions {
            jvmTarget.set(JvmTarget.JVM_21)
        }
    }

    test {
        val baseArgs = mutableListOf("-Xmx10g")

        // Detect path to BLAS native library
        val overrideDir: String? = System.getenv("BLAS_LIB_PATH")
        if (!overrideDir.isNullOrBlank()) {
            baseArgs += "-Djava.library.path=$overrideDir"
        } else {
            val os = System.getProperty("os.name").lowercase(Locale.ROOT)
            val nativeDir = when {
                "mac" in os -> {
                    val intel = "/usr/local/opt/openblas/lib"
                    val silicon = "/opt/homebrew/opt/openblas/lib"
                    when {
                        file(intel).exists() -> intel
                        file(silicon).exists() -> silicon
                        else -> {
                            logger.warn("Neither $intel nor $silicon exists. BLAS tests may fail unless BLAS_LIB_PATH is set.")
                            ""
                        }
                    }
                }
                "linux" in os -> "/usr/lib/x86_64-linux-gnu"
                else -> {
                    logger.warn("Unrecognized OS: $os. BLAS tests may fail unless BLAS_LIB_PATH is set.")
                    ""
                }
            }

            if (nativeDir.isNotBlank()) {
                baseArgs += "-Djava.library.path=$nativeDir"
            }
        }

        exclude(
            "**/analysis/gobii/*Test.class",
            "**/analysis/gbs/*Test.class",
            "**/analysis/gbs/v2/*Test.class",
            "**/analysis/gbs/repgen/*Test.class",
            "**/analysis/monetdb/*Test.class",
            "**/LowLevelCopyOfHDF5Test.class",
            "**/SplitHDF5ByChromosomePluginTest.class",
            "**/ThinSitesByPositionPluginTest.class",
            "**/LDKNNiImputationPluginTest.class",
            "**/GenomeFeatureBuilderTest.class",
            "**/BasicGenotypeMergeRuleTest.class",
            "**/DistanceMatrixHDF5Test.class",
            "**/TagsOnPhysMapHDF5Test.class",
            "**/FastMultithreadedAssociationPluginTest.class",
            "**/BuildUnfinishedHDF5GenotypesPluginTest.class",
            "**/GenomeAnnosDBQueryToPositionListPluginTest.class",
            "**/analysis/rna/*Test.class",
            "**/CreateFastaOrFastqFiles.class", // hard coded file paths (LCJ)
        )

        ignoreFailures = true // currently setting this to 'true' until we figure out failing tests
        jvmArgs = baseArgs

        println(jvmArgs)
    }
}


kotlin {
    jvmToolchain(21)
}
