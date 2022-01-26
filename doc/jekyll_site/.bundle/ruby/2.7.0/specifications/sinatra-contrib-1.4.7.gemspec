# -*- encoding: utf-8 -*-
# stub: sinatra-contrib 1.4.7 ruby lib

Gem::Specification.new do |s|
  s.name = "sinatra-contrib".freeze
  s.version = "1.4.7"

  s.required_rubygems_version = Gem::Requirement.new(">= 0".freeze) if s.respond_to? :required_rubygems_version=
  s.require_paths = ["lib".freeze]
  s.authors = ["Konstantin Haase".freeze, "Zachary Scott".freeze, "Gabriel Andretta".freeze, "Trevor Bramble".freeze, "Katrina Owen".freeze, "Ashley Williams".freeze, "Nicolas Sanguinetti".freeze, "Hrvoje \u0160imi\u0107".freeze, "Masahiro Fujiwara".freeze, "Rafael Magana".freeze, "Vipul A M".freeze, "ashley williams".freeze, "Jack Chu".freeze, "Sumeet Singh".freeze, "Ilya Shindyapin".freeze, "lest".freeze, "Jake Worth".freeze, "Kashyap".freeze, "Matt Lyon".freeze, "Matthew Conway".freeze, "Meck".freeze, "Michi Huber".freeze, "Nic Benders".freeze, "Patricio Mac Adden".freeze, "Reed Lipman".freeze, "Samy Dindane".freeze, "Sergey Nartimov".freeze, "Thibaut Sacreste".freeze, "Uchio KONDO".freeze, "Will Bailey".freeze, "mono".freeze, "Adrian Paca\u0142a".freeze, "undr".freeze, "Aish".freeze, "Alexey Chernenkov".freeze, "Andrew Crump".freeze, "Anton Davydov".freeze, "Bo Jeanes".freeze, "David Asabina".freeze, "Eliot Shepard".freeze, "Eric Marden".freeze, "Gray Manley".freeze, "Guillaume Bouteille".freeze, "Jamie Hodge".freeze, "Julie Ng".freeze, "Koichi Sasada".freeze, "Kyle Lacy".freeze, "Lars Vonk".freeze, "Martin Frost".freeze, "Mathieu Allaire".freeze]
  s.date = "2016-04-11"
  s.description = "Collection of useful Sinatra extensions".freeze
  s.email = ["konstantin.mailinglists@googlemail.com".freeze, "ohhgabriel@gmail.com".freeze, "inbox@trevorbramble.com".freeze, "e@zzak.io".freeze, "zachary@zacharyscott.net".freeze, "katrina.owen@gmail.com".freeze, "ashley@bocoup.com".freeze, "contacto@nicolassanguinetti.info".freeze, "shime.ferovac@gmail.com".freeze, "raf.magana@gmail.com".freeze, "m-fujiwara@axsh.net".freeze, "vipulnsward@gmail.com".freeze, "konstantin.haase@gmail.com".freeze, "jack@jackchu.com".freeze, "ashley666ashley@gmail.com".freeze, "ilya@shindyapin.com".freeze, "just.lest@gmail.com".freeze, "kashyap.kmbc@gmail.com".freeze, "ortuna@gmail.com".freeze, "tbramble@chef.io".freeze, "jworth@prevailhs.com".freeze, "mail@zzak.io".freeze, "nic@newrelic.com".freeze, "patriciomacadden@gmail.com".freeze, "rmlipman@gmail.com".freeze, "samy@dindane.com".freeze, "just.lest@gmail.com".freeze, "thibaut.sacreste@gmail.com".freeze, "udzura@udzura.jp".freeze, "will.bailey@gmail.com".freeze, "mono@mono0x.net".freeze, "altpacala@gmail.com".freeze, "undr@yandex.ru".freeze, "aisha.fenton@visfleet.com".freeze, "laise@pisem.net".freeze, "andrew.crump@ieee.org".freeze, "antondavydov.o@gmail.com".freeze, "me@bjeanes.com".freeze, "david@supr.nu".freeze, "eshepard@slower.net".freeze, "eric.marden@gmail.com".freeze, "g.manley@tukaiz.com".freeze, "duffman@via.ecp.fr".freeze, "jamiehodge@me.com".freeze, "uxjulie@gmail.com".freeze, "ko1@atdot.net".freeze, "kylewlacy@me.com".freeze, "lars.vonk@gmail.com".freeze, "blame@kth.se".freeze, "mathieuallaire@gmail.com".freeze, "matt@flowerpowered.com".freeze, "himself@mattonrails.com".freeze, "yesmeck@gmail.com".freeze, "michi.huber@gmail.com".freeze]
  s.homepage = "http://github.com/sinatra/sinatra-contrib".freeze
  s.licenses = ["MIT".freeze]
  s.rubygems_version = "3.2.27".freeze
  s.summary = "Collection of useful Sinatra extensions".freeze

  s.installed_by_version = "3.2.27" if s.respond_to? :installed_by_version

  if s.respond_to? :specification_version then
    s.specification_version = 4
  end

  if s.respond_to? :add_runtime_dependency then
    s.add_runtime_dependency(%q<sinatra>.freeze, ["~> 1.4.0"])
    s.add_runtime_dependency(%q<backports>.freeze, [">= 2.0"])
    s.add_runtime_dependency(%q<tilt>.freeze, [">= 1.3", "< 3"])
    s.add_runtime_dependency(%q<rack-test>.freeze, [">= 0"])
    s.add_runtime_dependency(%q<rack-protection>.freeze, [">= 0"])
    s.add_runtime_dependency(%q<multi_json>.freeze, [">= 0"])
    s.add_development_dependency(%q<rspec>.freeze, ["~> 2.3"])
    s.add_development_dependency(%q<haml>.freeze, [">= 0"])
    s.add_development_dependency(%q<erubis>.freeze, [">= 0"])
    s.add_development_dependency(%q<slim>.freeze, [">= 0"])
    s.add_development_dependency(%q<less>.freeze, [">= 0"])
    s.add_development_dependency(%q<sass>.freeze, [">= 0"])
    s.add_development_dependency(%q<builder>.freeze, [">= 0"])
    s.add_development_dependency(%q<liquid>.freeze, [">= 0"])
    s.add_development_dependency(%q<redcarpet>.freeze, [">= 0"])
    s.add_development_dependency(%q<RedCloth>.freeze, [">= 0"])
    s.add_development_dependency(%q<asciidoctor>.freeze, [">= 0"])
    s.add_development_dependency(%q<radius>.freeze, [">= 0"])
    s.add_development_dependency(%q<coffee-script>.freeze, [">= 0"])
    s.add_development_dependency(%q<nokogiri>.freeze, [">= 0"])
    s.add_development_dependency(%q<creole>.freeze, [">= 0"])
    s.add_development_dependency(%q<wikicloth>.freeze, [">= 0"])
    s.add_development_dependency(%q<markaby>.freeze, [">= 0"])
    s.add_development_dependency(%q<rake>.freeze, ["< 11"])
  else
    s.add_dependency(%q<sinatra>.freeze, ["~> 1.4.0"])
    s.add_dependency(%q<backports>.freeze, [">= 2.0"])
    s.add_dependency(%q<tilt>.freeze, [">= 1.3", "< 3"])
    s.add_dependency(%q<rack-test>.freeze, [">= 0"])
    s.add_dependency(%q<rack-protection>.freeze, [">= 0"])
    s.add_dependency(%q<multi_json>.freeze, [">= 0"])
    s.add_dependency(%q<rspec>.freeze, ["~> 2.3"])
    s.add_dependency(%q<haml>.freeze, [">= 0"])
    s.add_dependency(%q<erubis>.freeze, [">= 0"])
    s.add_dependency(%q<slim>.freeze, [">= 0"])
    s.add_dependency(%q<less>.freeze, [">= 0"])
    s.add_dependency(%q<sass>.freeze, [">= 0"])
    s.add_dependency(%q<builder>.freeze, [">= 0"])
    s.add_dependency(%q<liquid>.freeze, [">= 0"])
    s.add_dependency(%q<redcarpet>.freeze, [">= 0"])
    s.add_dependency(%q<RedCloth>.freeze, [">= 0"])
    s.add_dependency(%q<asciidoctor>.freeze, [">= 0"])
    s.add_dependency(%q<radius>.freeze, [">= 0"])
    s.add_dependency(%q<coffee-script>.freeze, [">= 0"])
    s.add_dependency(%q<nokogiri>.freeze, [">= 0"])
    s.add_dependency(%q<creole>.freeze, [">= 0"])
    s.add_dependency(%q<wikicloth>.freeze, [">= 0"])
    s.add_dependency(%q<markaby>.freeze, [">= 0"])
    s.add_dependency(%q<rake>.freeze, ["< 11"])
  end
end
