# frozen_string_literal: true

module Mopti
  # ScaledConjugateGradient is a class that implements multivariate optimization using scaled conjugate gradient method.
  #
  # @example
  #   # Seek the minimum value of the expression a*u**2 + b*u*v + c*v**2 + d*u + e*v + f for
  #   # given values of the parameters and an initial guess (u, v) = (0, 0).
  #   # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_cg.html#scipy.optimize.fmin_cg
  #   require 'numo/narray'
  #   require 'mopti'
  #
  #   args = [2, 3, 7, 8, 9, 10]
  #
  #   f = proc do |x, a, b, c, d, e, f|
  #     u = x[0]
  #     v = x[1]
  #     a * u**2 + b * u * v + c * v**2 + d * u + e * v + f
  #   end
  #
  #   g = proc do |x, a, b, c, d, e, f|
  #     u = x[0]
  #     v = x[1]
  #     gu = 2 * a * u + b * v + d
  #     gv = b * u + 2 * c * v + e
  #     Numo::DFloat[gu, gv]
  #   end
  #
  #   x0 = Numo::DFloat.zeros(2)
  #
  #   optimizer = Mopti::ScaledConjugateGradient.new(fnc: f, jcb: g, x_init: x0, args: args)
  #   result = optimizer.map { |params| params }.last
  #
  #   pp result
  #
  #   # {:x=>Numo::DFloat#shape=[2]
  #   # [-1.80847, -0.25533],
  #   #  :n_fev=>10,
  #   #  :n_jev=>18,
  #   #  :n_iter=>9,
  #   #  :fnc=>1.6170212789006122,
  #   #  :jcb=>1.8698188678645777e-07}
  #
  # *Reference*
  # 1. Moller, M F., "A Scaled Conjugate Gradient Algorithm for Fast Supervised Learning," Neural Networks, Vol. 6, pp. 525--533, 1993.
  class ScaledConjugateGradient
    include Enumerable

    # Create a new optimizer with scaled conjugate gradient.
    #
    # @param fnc [Method/Proc] Method for calculating the objective function to be minimized.
    # @param jcb [Method/Proc] Method for calculating the gradient vector.
    # @param args [Array/Hash] Arguments pass to the 'fnc' and 'jcb'.
    # @param x_init [Numo::NArray] Initial point.
    # @param max_iter [Integer] Maximum number of iterations.
    # @param xtol [Float] Tolerance for termination for the optimal vector norm.
    # @param ftol [Float] Tolerance for termination for the objective function value.
    # @param jtol [Float] Tolerance for termination for the gradient norm.
    def initialize(fnc:, jcb:, x_init:, args: nil, max_iter: 200, xtol: 1e-6, ftol: 1e-8, jtol: 1e-7)
      @fnc = fnc
      @jcb = jcb
      @x_init = x_init
      @args = args
      @max_iter = max_iter
      @ftol = ftol
      @jtol = jtol
      @xtol = xtol
    end

    # Iteration for optimization.
    #
    # @overload each(&block) -> Object
    # @yield [Hash] { x:, n_fev:, n_jev:, n_iter:, fnc:, jcb: }
    #   - x [Numo::NArray] Updated vector by optimization.
    #   - n_fev [Interger] Number of calls of the objective function.
    #   - n_jev [Integer] Number of calls of the jacobian.
    #   - n_iter [Integer] Number of iterations.
    #   - fnc [Float] Value of the objective function.
    #   - jcb [Numo::Narray] Values of the jacobian
    # @return [Enumerator] If block is not given, this method returns Enumerator.
    def each
      return to_enum(__method__) unless block_given?

      x = @x_init
      f_prev = func(x, @args)
      n_fev = 1
      f_curr = f_prev
      j_next = jacb(x, @args)
      n_jev = 1

      j_curr = j_next.dot(j_next)
      j_prev = j_next.dup
      d = -j_next
      success = true
      n_successes = 0
      beta = 1.0

      n_iter = 0

      while n_iter < @max_iter
        if success
          mu = d.dot(j_next)
          if mu >= 0.0
            d = -j_next
            mu = d.dot(j_next)
          end
          kappa = d.dot(d)
          break if kappa < 1e-16

          sigma = SIGMA_INIT / Math.sqrt(kappa)
          x_plus = x + (sigma * d)
          j_plus = jacb(x_plus, @args)
          n_jev += 1
          theta = d.dot(j_plus - j_next) / sigma
        end

        delta = theta + (beta * kappa)
        if delta <= 0
          delta = beta * kappa
          # TODO: Investigate the cause of the type error.
          # Cannot assign a value of type `::Complex` to a variable of type `::Float`
          # beta -= theta / kappa
          beta = (beta - (theta / kappa)).to_f
        end
        alpha = -mu / delta

        x_next = x + (alpha * d)
        f_next = func(x_next, @args)
        n_fev += 1

        delta = 2 * (f_next - f_prev) / (alpha * mu)
        if delta >= 0
          success = true
          n_successes += 1
          x = x_next
          f_curr = f_next
        else
          success = false
          f_curr = f_prev
        end

        n_iter += 1
        yield({ x: x, n_fev: n_fev, n_jev: n_jev, n_iter: n_iter, fnc: f_curr, jcb: j_curr })

        if success
          break if (f_next - f_prev).abs < @ftol
          break if (alpha * d).abs.max < @xtol

          f_prev = f_next

          j_prev = j_next
          j_next = jacb(x, @args)
          n_jev += 1

          j_curr = j_next.dot(j_next)
          break if j_curr <= @jtol
        end

        beta = [beta * 4, BETA_MAX].min if delta < 0.25
        beta = [beta / 4, BETA_MIN].max if delta > 0.75

        if n_successes == x.size
          d = -j_next
          beta = 1.0
          n_successes = 0
        elsif success
          gamma = (j_prev - j_next).dot(j_next) / mu
          d = -j_next + (gamma * d)
        end
      end
    end

    SIGMA_INIT = 1e-4
    BETA_MIN = 1e-15
    BETA_MAX = 1e+15

    private_constant :SIGMA_INIT, :BETA_MIN, :BETA_MAX

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

    def jacb(x, args)
      if args.is_a?(Hash)
        @jcb.call(x, **args)
      elsif args.is_a?(Array)
        @jcb.call(x, *args)
      elsif args.nil?
        @jcb.call(x)
      else
        @jcb.call(x, args)
      end
    end
  end
end
