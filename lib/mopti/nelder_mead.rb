# frozen_string_literal: true

module Mopti
  # NelderMead is a class that implements multivariate optimization using the Nelder-Mead simplex method.
  #
  # @example
  #   require 'numo/narray'
  #   require 'mopti'
  #
  #   args = Numo::DFloat[2, 3]
  #
  #   f = proc { |x, a| (8 * (x - a)**2).sum }
  #
  #   x0 = Numo::DFloat.zeros(2)
  #
  #   optimizer = Mopti::NelderMead.new(fnc: f, x_init: x0, args: args)
  #   result = optimizer.map { |params| params }.last
  #
  #   pp result
  #
  #   # {:x=>Numo::DFloat(view)#shape=[2]
  #   # [2, 3],
  #   #  :n_fev=>165,
  #   #  :n_iter=>84,
  #   #  :fnc=>5.694864987422661e-13}
  #
  # *Reference*
  # 1. Gao, F. and Han, L., "Implementing the Nelder-Mead simplex algorithm with adaptive parameters," Computational Optimization and Applications, vol. 51 (1), pp. 259--277, 2012.
  class NelderMead
    include Enumerable

    # Create a new optimizer with the Nelder-Mead simplex method.
    #
    # @param fnc [Method/Proc] Method for calculating the objective function to be minimized.
    # @param args [Array/Hash] Arguments pass to the 'fnc' and 'jcb'.
    # @param x_init [Numo::NArray] Initial point.
    # @param max_iter [Integer] Maximum number of iterations.
    #   If nil is given, max_iter sets to 200 * number of dimensions.
    # @param xtol [Float] Tolerance for termination for the optimal vector norm.
    # @param ftol [Float] Tolerance for termination for the objective function value.
    def initialize(fnc:, x_init:, args: nil, max_iter: nil, xtol: 1e-6, ftol: 1e-6)
      @fnc = fnc
      @args = args
      @x_init = x_init
      @max_iter = max_iter
      @xtol = xtol
      @ftol = ftol
    end

    # Iteration for optimization.
    #
    # @overload each(&block) -> Object
    # @yield [Hash] { x:, n_fev:, n_jev:, n_iter:, fnc:, jcb: }
    #   - x [Numo::NArray] Updated vector by optimization.
    #   - n_fev [Interger] Number of calls of the objective function.
    #   - n_iter [Integer] Number of iterations.
    #   - fnc [Float] Value of the objective function.
    # @return [Enumerator] If block is not given, this method returns Enumerator.
    def each
      return to_enum(__method__) unless block_given?

      x = @x_init.dup
      n = x.size
      max_iter = @max_iter || (200 * n)

      alpha = 1.0
      beta = n > 1 ? 1 + 2.fdiv(n) : 2.0
      gamma = n > 1 ? 0.75 - 1.fdiv(2 * n) : 0.5
      delta = n > 1 ? 1 - 1.fdiv(n) : 0.5

      sim = x.class.zeros(n + 1, n)
      sim[0, true] = x
      n.times do |k|
        y = x.dup
        y[k] = (y[k]).zero? ? ZERO_TAU : (1 + NON_ZERO_TAU) * y[k]
        sim[k + 1, true] = y
      end

      fsim = Numo::DFloat.zeros(n + 1)

      (n + 1).times { |k| fsim[k] = func(sim[k, true], @args) }
      n_fev = n + 1

      ind = fsim.sort_index
      fsim = fsim[ind].dup
      sim = sim[ind, true].dup

      n_iter = 0
      while n_iter < max_iter
        break if ((sim[1..-1, true] - sim[0, true]).abs.flatten.max <= @xtol) && ((fsim[0] - fsim[1..-1]).abs.max <= @ftol)

        xbar = sim[0...-1, true].sum(axis: 0) / n
        xr = xbar + (alpha * (xbar - sim[-1, true]))
        fr = func(xr, @args)
        n_fev += 1

        shrink = true
        if fr < fsim[0]
          xe = xbar + (beta * (xr - xbar))
          fe = func(xe, @args)
          n_fev += 1
          shrink = false
          if fe < fr
            sim[-1, true] = xe
            fsim[-1] = fe
          else
            sim[-1, true] = xr
            fsim[-1] = fr
          end
        elsif fr < fsim[-2]
          shrink = false
          sim[-1, true] = xr
          fsim[-1] = fr
        elsif fr < fsim[-1]
          xoc = xbar + (gamma * (xr - xbar))
          foc = func(xoc, @args)
          n_fev += 1
          if foc <= fr
            shrink = false
            sim[-1, true] = xoc
            fsim[-1] = foc
          end
        else
          xic = xbar - (gamma * (xr - xbar))
          fic = func(xic, @args)
          n_fev += 1
          if fic < fsim[-1]
            shrink = false
            sim[-1, true] = xic
            fsim[-1] = fic
          end
        end

        if shrink
          (1..n).to_a.each do |j|
            sim[j, true] = sim[0, true] + (delta * (sim[j, true] - sim[0, true]))
            fsim[j] = func(sim[j, true], @args)
            n_fev += 1
          end
        end

        ind = fsim.sort_index
        sim = sim[ind, true].dup
        fsim = fsim[ind].dup

        n_iter += 1

        yield({ x: sim[0, true], n_fev: n_fev, n_iter: n_iter, fnc: fsim.min })
      end
    end

    NON_ZERO_TAU = 0.05
    ZERO_TAU = 0.00025

    private_constant :NON_ZERO_TAU, :ZERO_TAU

    private

    def func(x, args)
      if args.is_a?(Hash)
        @fnc.call(x, **args)
      elsif args.is_a?(Array)
        @fnc.call(x, *args)
      elsif args.nil?
        @fnc.call(x)
      else
        @fnc.call(x, args)
      end
    end
  end
end
