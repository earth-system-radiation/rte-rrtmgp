# -*- encoding: utf-8 -*-
# stub: jekyll-github-metadata 2.0.2 ruby lib

Gem::Specification.new do |s|
  s.name = "jekyll-github-metadata".freeze
  s.version = "2.0.2"

  s.required_rubygems_version = Gem::Requirement.new(">= 0".freeze) if s.respond_to? :required_rubygems_version=
  s.require_paths = ["lib".freeze]
  s.authors = ["Parker Moore".freeze]
  s.date = "2016-06-22"
  s.email = ["parkrmoore@gmail.com".freeze]
  s.homepage = "https://github.com/parkr/github-metadata".freeze
  s.licenses = ["MIT".freeze]
  s.rubygems_version = "3.2.27".freeze
  s.summary = "The site.github namespace".freeze

  s.installed_by_version = "3.2.27" if s.respond_to? :installed_by_version

  if s.respond_to? :specification_version then
    s.specification_version = 4
  end

  if s.respond_to? :add_runtime_dependency then
    s.add_runtime_dependency(%q<octokit>.freeze, ["~> 4.0"])
    s.add_runtime_dependency(%q<jekyll>.freeze, ["~> 3.1"])
    s.add_development_dependency(%q<bundler>.freeze, ["~> 1.5"])
    s.add_development_dependency(%q<rake>.freeze, [">= 0"])
    s.add_development_dependency(%q<netrc>.freeze, [">= 0"])
    s.add_development_dependency(%q<rspec>.freeze, [">= 0"])
  else
    s.add_dependency(%q<octokit>.freeze, ["~> 4.0"])
    s.add_dependency(%q<jekyll>.freeze, ["~> 3.1"])
    s.add_dependency(%q<bundler>.freeze, ["~> 1.5"])
    s.add_dependency(%q<rake>.freeze, [">= 0"])
    s.add_dependency(%q<netrc>.freeze, [">= 0"])
    s.add_dependency(%q<rspec>.freeze, [">= 0"])
  end
end
