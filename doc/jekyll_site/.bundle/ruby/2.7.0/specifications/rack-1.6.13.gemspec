# -*- encoding: utf-8 -*-
# stub: rack 1.6.13 ruby lib

Gem::Specification.new do |s|
  s.name = "rack".freeze
  s.version = "1.6.13"

  s.required_rubygems_version = Gem::Requirement.new(">= 0".freeze) if s.respond_to? :required_rubygems_version=
  s.require_paths = ["lib".freeze]
  s.authors = ["Christian Neukirchen".freeze]
  s.date = "2020-02-08"
  s.description = "Rack provides a minimal, modular and adaptable interface for developing\nweb applications in Ruby.  By wrapping HTTP requests and responses in\nthe simplest way possible, it unifies and distills the API for web\nservers, web frameworks, and software in between (the so-called\nmiddleware) into a single method call.\n\nAlso see http://rack.github.io/.\n".freeze
  s.email = "chneukirchen@gmail.com".freeze
  s.executables = ["rackup".freeze]
  s.extra_rdoc_files = ["README.rdoc".freeze, "KNOWN-ISSUES".freeze, "HISTORY.md".freeze]
  s.files = ["HISTORY.md".freeze, "KNOWN-ISSUES".freeze, "README.rdoc".freeze, "bin/rackup".freeze]
  s.homepage = "http://rack.github.io/".freeze
  s.licenses = ["MIT".freeze]
  s.rubygems_version = "3.2.27".freeze
  s.summary = "a modular Ruby webserver interface".freeze

  s.installed_by_version = "3.2.27" if s.respond_to? :installed_by_version

  if s.respond_to? :specification_version then
    s.specification_version = 4
  end

  if s.respond_to? :add_runtime_dependency then
    s.add_development_dependency(%q<bacon>.freeze, [">= 0"])
    s.add_development_dependency(%q<rake>.freeze, [">= 0"])
  else
    s.add_dependency(%q<bacon>.freeze, [">= 0"])
    s.add_dependency(%q<rake>.freeze, [">= 0"])
  end
end
