# frozen_string_literal: true

RSpec.describe Mopti do
  describe '#minimize' do
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

    let(:jcb) do
      proc do |x, f, k|
        log_pdot = f.dot(x)
        log_z = Numo::NMath.log(Numo::NMath.exp(log_pdot).sum)
        p = Numo::NMath.exp(log_pdot - log_z)
        f.transpose.dot(p) - k
      end
    end

    let(:res) { described_class.minimize(algorithm: algorithm, **args) }
    let(:err) { (fnc.call(res[:x], f, k) - fnc.call(Numo::DFloat[0, -0.524869316, 0.487525860], f, k)).abs }

    shared_examples 'perform optimization' do
      it 'obtains optimized vector that minimizes the value of the objective function' do
        expect(err).to be < 1.0e-6
      end
    end

    context 'when optimizer is scaled conjugate gradient' do
      let(:algorithm) { 'SCG' }
      let(:args) { { fnc: fnc, jcb: jcb, x_init: x, args: [f, k] } }

      it_behaves_like 'perform optimization'
    end

    context 'when optimizer is Nelder-Mead simplex' do
      let(:algorithm) { 'Nelder-Mead' }
      let(:args) { { fnc: fnc, x_init: x, args: [f, k] } }

      it_behaves_like 'perform optimization'
    end

    context 'when non-existent optimizer is given' do
      it 'raises ArgumentError' do
        expect { described_class.minimize(algorithm: 'foobar', fnc: fnc) }.to raise_error(ArgumentError)
      end
    end
  end
end
