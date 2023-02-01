@file:JvmName("TasselVersions")

package net.maizegenetics.tassel

import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.util.Utils
import org.apache.log4j.Logger
import org.w3c.dom.Node
import javax.xml.parsers.DocumentBuilderFactory

object TasselVersions {

    private val myLogger = Logger.getLogger(LibraryInfo::class.java)

    private val infoMap = mutableMapOf<String, LibraryInfo>()

    data class LibraryInfo(val name: String, val version: String, val date: String, val citation: String?)

    init {
        for (current in Utils.getFullyQualifiedResourceNames("tassel_info.xml")) {
            getThirdPartyLibraryInfo(current)
        }
    }

    fun tasselName() = ""

    fun tasselVersion() = TASSELMainFrame.version

    fun tasselVersionDate() = TASSELMainFrame.versionDate

    fun tasselCitation() = AbstractPlugin.DEFAULT_CITATION

    fun phgName() = infoMap["PHG"]!!.name

    fun phgVersion() = infoMap["PHG"]!!.version

    fun phgVersionDate() = infoMap["PHG"]!!.date

    fun phgCitation() = infoMap["PHG"]!!.citation

    fun libraryInfo(library: String): LibraryInfo? = infoMap[library]

    fun libraryInfos() = infoMap.entries

    private fun getThirdPartyLibraryInfo(filename: String) {
        try {
            LibraryInfo::class.java.getResourceAsStream(filename).use { input ->
                val dbFactory = DocumentBuilderFactory.newInstance()
                val dBuilder = dbFactory.newDocumentBuilder()
                val doc = dBuilder.parse(input)
                doc.documentElement.normalize()
                val rootElement = doc.documentElement
                require(
                    rootElement.nodeName.equals(
                        "Tassel",
                        ignoreCase = true
                    )
                ) { """LibraryInfo: getThirdPartyLibraryInfo: Root Node must be Tassel: ${rootElement.nodeName}""" }
                val children = rootElement.childNodes
                for (i in 0 until children.length) {
                    val current = children.item(i)
                    addLibInfo(current)
                }
            }
        } catch (e: Exception) {
            myLogger.debug(e.message, e)
            throw IllegalStateException(
                """
                LibraryInfo: getThirdPartyLibraryInfo: Problem reading XML file: $filename
                ${e.message}
                """.trimIndent()
            )
        }
    }

    private fun addLibInfo(rootElement: Node) {

        try {
            if (!rootElement.nodeName.trim { it <= ' ' }.equals("TasselLibrary", ignoreCase = true)) {
                return
            }
            val children = rootElement.childNodes
            val libraryInfo = mutableMapOf<String, String>()
            for (i in 0 until children.length) {
                val current = children.item(i)
                val elementName = current.nodeName.trim()
                if (current.nodeType == Node.ELEMENT_NODE) {
                    libraryInfo[elementName] = current.textContent.trim()
                }
            }
            infoMap[libraryInfo["key"]!!] =
                LibraryInfo(
                    libraryInfo["name"]!!,
                    libraryInfo["version"]!!,
                    libraryInfo["date"]!!,
                    libraryInfo["citation"]!!
                )
        } catch (e: java.lang.Exception) {
            myLogger.debug(e.message, e)
        }

    }

}