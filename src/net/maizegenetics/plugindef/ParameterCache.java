package net.maizegenetics.plugindef;

import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.util.Optional;
import java.util.Properties;

/**
 * This class is for storing a cache of parameters retrieved from a {@link java.util.Properties}
 * file that is available to all plugins on the command line.
 * This is used when run_pipeline.pl -configParameters config.txt... is used.
 *
 * @author Terry Casstevens
 * Created January 05, 2018
 */
public class ParameterCache {

    private static final Logger myLogger = Logger.getLogger(ParameterCache.class);

    private static Properties CACHE = null;

    private ParameterCache() {
    }

    /**
     * Loads the parameter cache with the values in the given {@link java.util.Properties} file.
     *
     * @param filename file name
     */
    public static void load(String filename) {

        if (CACHE != null) {
            throw new IllegalStateException("ParameterCache: load: cache already loaded.");
        }

        CACHE = new Properties();
        try (BufferedReader reader = Utils.getBufferedReader(filename)) {
            CACHE.load(reader);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("ParameterCache: load: problem reading properties file: " + filename);
        }

        for (String key : CACHE.stringPropertyNames()) {
            myLogger.info("ParameterCache: key: " + key + " value: " + CACHE.getProperty(key));
        }

    }

    /**
     * Returns the value if any for the given plugin and parameter.  Value for PluginClassName.parameter will
     * be returned if it exists.  If not, value for parameter will be returned.  Otherwise, an empty optional.
     *
     * @param plugin plugin
     * @param parameter parameter
     *
     * @return value if exists
     */
    public static Optional<String> value(Plugin plugin, String parameter) {

        if (CACHE == null) {
            return Optional.empty();
        }

        String value = CACHE.getProperty(Utils.getBasename(plugin.getClass().getName()) + "." + parameter);
        if (value != null) {
            return Optional.of(value);
        }

        value = CACHE.getProperty(parameter);
        return Optional.ofNullable(value);

    }

    /**
     * Returns true if this cache has been loaded with values.
     *
     * @return true if this cache has been loaded with values.
     */
    public static boolean hasValues() {
        return CACHE != null;
    }

}
