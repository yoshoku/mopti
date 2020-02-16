# frozen_string_literal: true

require 'spec_helper'

RSpec.describe Mopti::NelderMead do
  let(:x) { Numo::DFloat.zeros(3) }
  let(:f) { Numo::DFloat[[1, 1, 1], [1, 1, 0], [1, 0, 1], [1, 0, 0], [1, 0, 0]] }
  let(:k) { Numo::DFloat[1.0, 0.3, 0.5] }

  let(:fnc) do
    proc do |x, f, k|
      log_pdot = f.dot(x)
      log_z = Numo::NMath.log(Numo::NMath.exp(log_pdot).sum)
      (log_z - k.dot(x)).to_f
    end
  end

  let(:optimizer) { described_class.new(fnc: fnc, x_init: x, args: [f, k]) }
  let(:params) { optimizer.map { |p| p } }
  let(:res) { params.last }
  let(:err) { (fnc.call(res[:x], f, k) - fnc.call(Numo::DFloat[0, -0.524869316, 0.487525860], f, k)).abs }

  it { expect(params.size).to be_positive }
  it { expect(res.keys).to match(%i[x n_fev n_iter fnc]) }
  it { expect(err).to be < 1.0e-6 }
end
