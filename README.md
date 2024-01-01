# Mopti

[![Build Status](https://github.com/yoshoku/mopti/actions/workflows/build.yml/badge.svg)](https://github.com/yoshoku/mopti/actions/workflows/build.yml)
[![Gem Version](https://badge.fury.io/rb/mopti.svg)](https://badge.fury.io/rb/mopti)
[![BSD 3-Clause License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://github.com/yoshoku/mopti/blob/main/LICENSE.txt)
[![Documentation](https://img.shields.io/badge/api-reference-blue.svg)](https://yoshoku.github.io/mopti/doc/)

Mopti is a multivariate optimization library in Ruby.
Mopti supports Nelder-Mead simplex method and Scaled Conjugate Gradient.

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'mopti'
```

And then execute:

    $ bundle install

Or install it yourself as:

    $ gem install mopti

## Documentation

- [Mopti API Documentation](https://yoshoku.github.io/mopti/doc/)

## Usage

Example 1. Linear Regression with Nelder-Mead simplex method

```ruby
require 'numo/narray'
require 'numo/gnuplot'
require 'mopti'

# Define objective function.
obj_fnc = proc do |w, x, y|
  n_samples = x.shape[0]
  ((y - x.dot(w))**2).sum.fdiv(n_samples)
end

# Explanatory variables.
x = Numo::DFloat[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Dependent variables.
y = 3 * x + 2

# Extend variable for intercept.
xb = Numo::DFloat.vstack([x, [1] * x.size]).transpose.dup

# Optimize parameter vectors.
res = Mopti::minimize(algorithm: 'Nelder-Mead', fnc: obj_fnc, x_init: Numo::DFloat.zeros(2), args: [xb, y])

# Output result.
pp res

a, b = res[:x].to_a
Numo.gnuplot do
  set(terminal: 'png')
  set(output: 'example1.png')
  plot [[0, 10], [a * x[0] + b, a * x[-1] + b], w: :lines, lw: 1, t: 'parameters'],
       [x, y, pt: 6, ps: 2, t: 'data']
end
```

```
$ brew install gnuplot
$ gem install mopti numo-narray numo-gnuplot
$ ruby example1.rb
{:x=>Numo::DFloat(view)#shape=[2]
[3, 2],
 :n_fev=>177,
 :n_iter=>91,
 :fnc=>2.290874014308807e-13}
```

![example1](https://user-images.githubusercontent.com/5562409/74737586-79c9db00-5298-11ea-9063-90e4655f878a.png)


## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `rake spec` to run the tests. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/yoshoku/mopti. This project is intended to be a safe, welcoming space for collaboration, and contributors are expected to adhere to the [code of conduct](https://github.com/yoshoku/mopti/blob/main/CODE_OF_CONDUCT.md).


## Code of Conduct

Everyone interacting in the Mopti project's codebases, issue trackers, chat rooms and mailing lists is expected to follow the [code of conduct](https://github.com/yoshoku/mopti/blob/main/CODE_OF_CONDUCT.md).
