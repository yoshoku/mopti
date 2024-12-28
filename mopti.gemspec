# frozen_string_literal: true

require_relative 'lib/mopti/version'

Gem::Specification.new do |spec|
  spec.name          = 'mopti'
  spec.version       = Mopti::VERSION
  spec.authors       = ['yoshoku']
  spec.email         = ['yoshoku@outlook.com']

  spec.summary       = 'Multivariate Optimization Library in Ruby.'
  spec.description   = <<~MSG
    Multivariate Optimization Library in Ruby.
    Mopti supports Nelder-Mead simplex method and Scaled Conjugate Gradient.
  MSG
  spec.homepage      = 'https://github.com/yoshoku/mopti'
  spec.license       = 'BSD-3-Clause'

  spec.metadata['homepage_uri'] = spec.homepage
  spec.metadata['source_code_uri'] = 'https://github.com/yoshoku/mopti'
  spec.metadata['changelog_uri'] = 'https://github.com/yoshoku/mopti/blob/master/CHANGELOG.md'
  spec.metadata['documentation_url'] = 'https://yoshoku.github.io/mopti/doc/'
  spec.metadata['rubygems_mfa_required'] = 'true'

  # Specify which files should be added to the gem when it is released.
  # The `git ls-files -z` loads the files in the RubyGem that have been added into git.
  spec.files = Dir.chdir(File.expand_path(__dir__)) do
    `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features|sig-deps)/}) }
  end
  spec.bindir        = 'exe'
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ['lib']

  spec.add_dependency 'numo-narray', '>= 0.9.1'
end
