# -*- encoding: utf-8 -*-
# stub: net-dns 0.9.0 ruby lib

Gem::Specification.new do |s|
  s.name = "net-dns".freeze
  s.version = "0.9.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0".freeze) if s.respond_to? :required_rubygems_version=
  s.require_paths = ["lib".freeze]
  s.authors = ["Marco Ceresa".freeze, "Simone Carletti".freeze]
  s.date = "2018-10-29"
  s.description = "Net::DNS is a pure Ruby DNS library, with a clean OO interface and an extensible API.".freeze
  s.email = ["ceresa@gmail.com".freeze, "weppos@weppos.net".freeze]
  s.extra_rdoc_files = ["LICENSE.txt".freeze]
  s.files = ["LICENSE.txt".freeze]
  s.homepage = "http://github.com/bluemonk/net-dns".freeze
  s.required_ruby_version = Gem::Requirement.new(">= 2.1".freeze)
  s.rubygems_version = "3.2.27".freeze
  s.summary = "Pure Ruby DNS library.".freeze

  s.installed_by_version = "3.2.27" if s.respond_to? :installed_by_version

  if s.respond_to? :specification_version then
    s.specification_version = 4
  end

  if s.respond_to? :add_runtime_dependency then
    s.add_development_dependency(%q<mocha>.freeze, [">= 0"])
    s.add_development_dependency(%q<rake>.freeze, [">= 0"])
    s.add_development_dependency(%q<yard>.freeze, [">= 0"])
  else
    s.add_dependency(%q<mocha>.freeze, [">= 0"])
    s.add_dependency(%q<rake>.freeze, [">= 0"])
    s.add_dependency(%q<yard>.freeze, [">= 0"])
  end
end
