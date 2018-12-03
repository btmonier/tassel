/*
 *  GeneratePluginCode
 */
package net.maizegenetics.plugindef;

import com.google.common.base.CaseFormat;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.awt.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;

/**
 *
 * @author Terry Casstevens
 */
public class GenerateRCode {

    private static final Logger myLogger = Logger.getLogger(GenerateRCode.class);

    private GenerateRCode() {

    }

    public static void main(String[] args) {
        generate(FilterSiteBuilderPlugin.class, "genotypeTable");
        //FilterSiteBuilderPlugin filterPlugin = new FilterSiteBuilderPlugin();
        //filterPlugin.setParameter("siteMinCount","0");
        //filterPlugin.runPlugin()
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

    private static void generate(AbstractPlugin plugin,String inputObject) {
        String clazz = Utils.getBasename(plugin.getClass().getName());
        StringBuilder sb=new StringBuilder("rTASSEL::"+stringToCamelCase(clazz));
        sb.append(" <- function(");
        sb.append(inputObject+",\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            //String guiNameAsCamelCase = stringToCamelCase(current.guiName());
            String methodName = removeMyFromString(field.getName());
            if(!(current.defaultValue() instanceof Number)) continue;
            sb.append("            "+methodName+"="+current.defaultValue()+",\n");

        }
        sb.deleteCharAt(sb.lastIndexOf(","));//remove last comma

        sb.append(") {\n");
        sb.append("    plugin <- rJava::.jnew(\""+plugin.getClass().getCanonicalName()+"\")\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            //String guiNameAsCamelCase = stringToCamelCase(current.guiName());
            if(!(current.defaultValue() instanceof Number)) continue;
            String methodName = removeMyFromString(field.getName());
            sb.append("    plugin$setParameter(\""+methodName+"\",toString("+methodName+"))\n");
        }
        sb.append("    plugin$runPlugin("+inputObject+")\n");
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
