/*
 *  GenerateRCode
 */
package net.maizegenetics.plugindef;

import com.google.common.base.CaseFormat;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.awt.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.Date;

/**
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class GenerateRCode {

    private static final Logger myLogger = Logger.getLogger(GenerateRCode.class);

    private GenerateRCode() {

    }

    public static void main(String[] args) {
        printHeader();
        generate(FilterSiteBuilderPlugin.class, "genotypeTable");
        generate(FilterTaxaBuilderPlugin.class, "genotypeTable");
        generate(KinshipPlugin.class, "genotypeTable");
    }

    private static void printHeader() {

        System.out.println("#!/usr/bin/env Rscript");

        System.out.println("\n#--------------------------------------------------------------------");
        System.out.println("# Script Name:   TasselPluginWrappers.R");
        System.out.println("# Description:   Generated R interface to TASSEL 5");
        System.out.println("# Author:        Brandon Monier, Ed Buckler, Terry Casstevens");
        System.out.print("# Created:       ");
        System.out.println(new Date());
        System.out.println("#--------------------------------------------------------------------");

        System.out.println("# Preamble\n");
        System.out.println("source(\"R/AllClasses.R\")");

        System.out.println("\n## Load packages");
        System.out.println("if (!requireNamespace(\"BiocManager\")) {");
        System.out.println("    install.packages(\"BiocManager\")");
        System.out.println("}");

        System.out.println("\npackages <- c(");
        System.out.println("\"rJava\"");
        System.out.println(")");
        System.out.println("BiocManager::install(packages)");
        System.out.println("library(rJava)");

        System.out.println("\n## Init JVM");
        System.out.println("rJava::.jinit()");

        System.out.println("\n## Add TASSEL 5 class path");
        System.out.println("rJava::.jaddClassPath(\"/tassel-5-standalone/lib\")");
        System.out.println("rJava::.jaddClassPath(\"/tassel-5-standalone/sTASSEL.jar\")\n");

    }

    public static void generate(Class currentMatch, String inputObject) {
        try {
            Constructor constructor = currentMatch.getConstructor(Frame.class);
            generate((AbstractPlugin) constructor.newInstance((Frame) null), inputObject);
        } catch (Exception ex) {
            try {
                Constructor constructor = currentMatch.getConstructor(Frame.class, boolean.class);
                generate((AbstractPlugin) constructor.newInstance(null, false), inputObject);
            } catch (NoSuchMethodException nsme) {
                myLogger.warn("Self-describing Plugins should implement this constructor: " + currentMatch.getClass().getName());
                myLogger.warn("public Plugin(Frame parentFrame, boolean isInteractive) {");
                myLogger.warn("   super(parentFrame, isInteractive);");
                myLogger.warn("}");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static void generate(AbstractPlugin plugin, String inputObject) {

        String clazz = Utils.getBasename(plugin.getClass().getName());
        //StringBuilder sb = new StringBuilder("rTASSEL::" + stringToCamelCase(clazz));
        StringBuilder sb = new StringBuilder(stringToCamelCase(clazz));
        sb.append(" <- function(");
        sb.append(inputObject + ",\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            String methodName = removeMyFromString(field.getName());
            if ((current.defaultValue() instanceof Number)) {
                sb.append("            " + methodName + "=" + current.defaultValue() + ",\n");
            } else if ((current.defaultValue() instanceof Enum)) {
                sb.append("            " + methodName + "=" + current.defaultValue() + ",\n");
            }

        }
        sb.deleteCharAt(sb.lastIndexOf(",")); //remove last comma

        sb.append(") {\n");
        sb.append("    plugin <- rJava::.jnew(\"" + plugin.getClass().getCanonicalName() + "\")\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            if ((current.defaultValue() instanceof Number)) {
                String methodName = removeMyFromString(field.getName());
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            } else if ((current.defaultValue() instanceof Enum)) {
                String methodName = removeMyFromString(field.getName());
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            }
        }
        sb.append("    plugin$runPlugin(" + inputObject + ")\n");
        sb.append("}\n");
        System.out.println(sb.toString());

    }

    private static final int DEFAULT_DESCRIPTION_LINE_LENGTH = 50;

    private static String createDescription(String description) {
        int count = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("     * ");
        for (int i = 0, n = description.length(); i < n; i++) {
            count++;
            if (description.charAt(i) == '\n') {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else if ((count > DEFAULT_DESCRIPTION_LINE_LENGTH) && (description.charAt(i) == ' ')) {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else {
                builder.append(description.charAt(i));
            }
        }
        return builder.toString();
    }

    private static String stringToCamelCase(String str) {
        StringBuilder builder = new StringBuilder();
        builder.append(Character.toLowerCase(str.charAt(0)));
        boolean makeUpper = false;
        for (int i = 1; i < str.length(); i++) {
            char current = str.charAt(i);
            if (current == ' ') {
                makeUpper = true;
            } else if (makeUpper) {
                builder.append(Character.toUpperCase(current));
                makeUpper = false;
            } else {
                builder.append(current);
            }
        }
        return builder.toString();
    }

    private static String removeMyFromString(String str) {
        String lower = str.toLowerCase();
        if (lower.startsWith("my")) {
            str = str.substring(2);
        }
        return CaseFormat.UPPER_CAMEL.to(CaseFormat.LOWER_CAMEL, str);
    }

}
