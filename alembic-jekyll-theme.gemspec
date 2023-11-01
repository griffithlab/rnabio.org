# coding: utf-8

Gem::Specification.new do |spec|
  spec.name          = "alembic-jekyll-theme"
  spec.version       = "1.5.2"
  spec.authors       = ["David Darnes"]
  spec.email         = ["me@daviddarnes.com"]

  spec.summary       = %q{A Jekyll boilerplate theme designed to be a starting point for any Jekyll website.}
  spec.description   = "A Jekyll boilerplate theme designed to be a starting point for any Jekyll website. Rather than starting from scratch, this boilerplate is designed to get the ball rolling immediately."
  spec.homepage      = "https://alembic.darn.es"
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0").select { |f| f.match(%r{^(assets|_layouts|_includes|_sass|LICENSE|README)}i) }

  spec.add_runtime_dependency "jekyll", "~> 4.1"
  spec.add_runtime_dependency "jekyll-sitemap", "~> 1.4.0"
  spec.add_runtime_dependency "jekyll-mentions", "~> 1.6.0"
  spec.add_runtime_dependency "jekyll-paginate", "~> 1.1.0"
  spec.add_runtime_dependency "jekyll-seo-tag", "~> 2.7.1"
  spec.add_runtime_dependency "jekyll-redirect-from", "~> 0.16"
  spec.add_runtime_dependency "jekyll-feed", "~> 0.15"
  spec.add_runtime_dependency "jemoji", "~> 0.12"
end
