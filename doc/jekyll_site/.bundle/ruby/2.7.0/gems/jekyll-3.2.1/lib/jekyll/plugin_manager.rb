module Jekyll
  class PluginManager
    attr_reader :site

    # Create an instance of this class.
    #
    # site - the instance of Jekyll::Site we're concerned with
    #
    # Returns nothing
    def initialize(site)
      @site = site
    end

    # Require all the plugins which are allowed.
    #
    # Returns nothing
    def conscientious_require
      require_plugin_files
      require_gems
      deprecation_checks
    end

    # Require each of the gem plugins specified.
    #
    # Returns nothing.
    def require_gems
      Jekyll::External.require_with_graceful_fail(
        site.gems.select { |gem| plugin_allowed?(gem) }
      )
    end

    def self.require_from_bundler
      if !ENV["JEKYLL_NO_BUNDLER_REQUIRE"] && File.file?("Gemfile")
        require "bundler"

        Bundler.setup
        required_gems = Bundler.require(:jekyll_plugins)
        message = "Required #{required_gems.map(&:name).join(", ")}"
        Jekyll.logger.debug("PluginManager:", message)
        ENV["JEKYLL_NO_BUNDLER_REQUIRE"] = "true"

        true
      else
        false
      end
    end

    # Check whether a gem plugin is allowed to be used during this build.
    #
    # gem_name - the name of the gem
    #
    # Returns true if the gem name is in the whitelist or if the site is not
    #   in safe mode.
    def plugin_allowed?(gem_name)
      !site.safe || whitelist.include?(gem_name)
    end

    # Build an array of allowed plugin gem names.
    #
    # Returns an array of strings, each string being the name of a gem name
    #   that is allowed to be used.
    def whitelist
      @whitelist ||= Array[site.config["whitelist"]].flatten
    end

    # Require all .rb files if safe mode is off
    #
    # Returns nothing.
    def require_plugin_files
      unless site.safe
        plugins_path.each do |plugin_search_path|
          plugin_files = Utils.safe_glob(plugin_search_path, File.join("**", "*.rb"))
          Jekyll::External.require_with_graceful_fail(plugin_files)
        end
      end
    end

    # Public: Setup the plugin search path
    #
    # Returns an Array of plugin search paths
    def plugins_path
      if site.config["plugins_dir"].eql? Jekyll::Configuration::DEFAULTS["plugins_dir"]
        [site.in_source_dir(site.config["plugins_dir"])]
      else
        Array(site.config["plugins_dir"]).map { |d| File.expand_path(d) }
      end
    end

    def deprecation_checks
      pagination_included = (site.config["gems"] || []).include?("jekyll-paginate") ||
        defined?(Jekyll::Paginate)
      if site.config["paginate"] && !pagination_included
        Jekyll::Deprecator.deprecation_message "You appear to have pagination " \
          "turned on, but you haven't included the `jekyll-paginate` gem. " \
          "Ensure you have `gems: [jekyll-paginate]` in your configuration file."
      end
    end
  end
end
