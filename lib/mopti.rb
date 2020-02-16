# frozen_string_literal: true

require 'numo/narray'

require 'mopti/version'
require 'mopti/nelder_mead'
require 'mopti/scaled_conjugate_gradient'

module Mopti
  module_function

  # Perform minization of the objective function.
  #
  # @param algorithm [String] Type of optimizer.
  #   - 'SCG': ScaledConjugateGradient
  #   - 'Nelder-Mead': NelderMead
  # @return [Hash] Result of optimization.
  def minimize(algorithm:, **args)
    optimizer = case algorithm
                when 'SCG'
                  ScaledConjugateGradient.new(**args)
                when 'Nelder-Mead'
                  NelderMead.new(**args)
                else
                  raise ArgumentError, 'A non-existent algorithm is specified'
                end
    res = nil
    optimizer.each { |params| res = params }
    res
  end
end
