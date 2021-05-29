# frozen_string_literal: true

source "https://rubygems.org"
gemspec

gem "github-pages", group: :jekyll_plugins

gem 'wdm', '>= 0.1.0'

gem "jekyll", ENV["JEKYLL_VERSION"] if ENV["JEKYLL_VERSION"]

install_if -> { Gem.win_platform? } do
  gem "tzinfo", "~> 1.2"
  gem "tzinfo-data"
end