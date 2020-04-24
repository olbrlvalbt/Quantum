classdef Quantum
    properties (Constant)
        I = [1 0; 0 1];
        X = [0 1; 1 0];
        Y = [0 -i; i 0];
        Z = [1 0; 0 -1];
        
        H = (1/sqrt(2)) * [1 1; 1 -1];

        u0 = [1; 0];
        u1 = [0; 1];
        
        up = Quantum.H * Quantum.u0;
        um = Quantum.H * Quantum.u1;

        u00 = kron(Quantum.u0, Quantum.u0);
        u01 = kron(Quantum.u0, Quantum.u1);
        u10 = kron(Quantum.u1, Quantum.u0);
        u11 = kron(Quantum.u1, Quantum.u1);
        
        P0 = Quantum.u0 * Quantum.u0';
        P1 = Quantum.u1 * Quantum.u1';
        Pp = Quantum.up * Quantum.up';
        Pm = Quantum.um * Quantum.um';

        u000 = kron(Quantum.u0, Quantum.u00);
        u001 = kron(Quantum.u0, Quantum.u01);
        u010 = kron(Quantum.u0, Quantum.u10);
        u011 = kron(Quantum.u0, Quantum.u11);
        u100 = kron(Quantum.u1, Quantum.u00);
        u101 = kron(Quantum.u1, Quantum.u01);
        u110 = kron(Quantum.u1, Quantum.u10);
        u111 = kron(Quantum.u1, Quantum.u11);

        b00 = (Quantum.u00 + Quantum.u11) / sqrt(2);
        b01 = (Quantum.u01 + Quantum.u10) / sqrt(2);
        b10 = (Quantum.u00 - Quantum.u11) / sqrt(2);
        b11 = (Quantum.u01 - Quantum.u10) / sqrt(2);
        
%         CN = kron(Quantum.P0, Quantum.I) + kron(Quantum.P1, Quantum.X);
    end
    methods
        function Pauli = getPauliMatrices(Quantum)
            Pauli = Quantum.I;
            Pauli(:, :, 2) = Quantum.X;
            Pauli(:, :, 3) = Quantum.Y;
            Pauli(:, :, 4) = Quantum.Z;
        end
        
        % Pure states
        
        function psi = createState1(Quantum, alpha, beta)
            psi = Quantum.normalizeState(alpha * Quantum.u0 + beta * Quantum.u1);
        end
        
        function psi = createState2(Quantum, alpha, beta, gamma, delta)
            psi = normalizeState(alpha * Quantum.u00 + beta * Quantum.u01 + gamma * Quantum.u10 + delta * Quantum.u11);
        end
        function CU = createControlledGate2(Quantum, U)
            CU = kron(Quantum.P0, Quantum.I) + kron(Quantum.P1, U);
        end
        
        function newPsi = normalizeState(Quantum, psi)
            newPsi = psi / sqrt(psi' * psi);
        end
        function rho = getDensityOperator(Quantum, psi)
            rho = psi * psi';
        end
        function [newPsi, pr] = measureState(Quantum, M, psi)
            newPsi = normalizeState(M * psi);
            pr = psi' * (M' * M) * psi;
        end
        function [newRho, pr] = measureDensityOperator(Quantum, M, rho)
            newRho = (M * rho * M') / trace((M' * M) * rho);
            pr = trace((M' * M) * psi);
        end
        function U = getBasisChangeMatrix(Quantum, originalBasis, newBasis)
%             if ~diff(size(originalBasis)) || ~diff(size(newBasis)) || ~diff(size(originalBasis), size(newBasis))
%                 throw(MException('Incorrect basis dimentions'));
%             end
            U = newBasis' * originalBasis;
        end
        function newPsi = changeBasisForState(Quantum, originalBasis, newBasis, psi)
            U = Quantum.getBasisChangeMatrix(originalBasis, newBasis);
            newPsi = U * psi
        end
        function newA = changeBasisForOperator(Quantum, originalBasis, newBasis, A)
            U = Quantum.getBasisChangeMatrix(originalBasis, newBasis);
            newA = U * A * U';
        end
        function expA = getExpectedValue(Quantum, A, psi)
            expA = psi' * A * psi; 
        end
        function uncertainty = getUncertainty(Quantum, A, psi)
            uncertainty = sqrt(Quantum.getExpectedValue(A * A, psi) - Quantum.getExpectedValue(A, psi)^2)
        end
        function commutator = commute(Quantum, A, B)
            commutator = A * B - B * A;
        end
        function deg = getDegeneration(Quantum, A, z)
            eigA = eig(A);
            deg = sum(eigA(:) == z);
        end
        function [Q, U, P] = getPolarDecomposition(Quantum, A)
            Q = sqrt(A * A');
            P = sqrt(A' * A);
            U = A / P;
        end
        
        function res = isPositiveDefinite(Quantum, A)
            res = isequal(A, A') && all(eig(A) >= 0);
        end
        function res = isPositiveSemidefinite(Quantum, A)
            res = isequal(A, A') && all(eig(A) > 0)
        end
        function res = doCommute(Quantum, A, B)
            commutator = Quantum.commute(A, B);
            res = all(commutator(:) == 0);
        end
        function res = isDegenerated(Quantum, A, z)
            res = Quantum.getDegeneration(A, z) - 1;
            if res > 0
                res = 1;
            end
        end
        
        
        % Mixed States
        
        function psi = createMixedState(Quantum, states, numberOfSystems)
            psi = MixedState(states, numberOfSystems);
            for i = 1:size(psi.States)
                psi.States(:,i) = Quantum.normalizeState(psi.States(:,i));
            end
        end
        function rho = getMixedDensityOperator(Quantum, psi)
            rho = 0;
            for i = 1:size(psi.States)
                col = psi.States(:,i);
                rho = rho + (col * col') * psi.getProbability(i);
            end
        end
        
        function res = isMixedDensityOperator(Quantum, rho)
            res = trace(rho * rho) < 1;
        end
        function res = isCompletelyMixedDensityOperator(Quantum, rho)
            [n, ] = size(rho)
            res = trace(rho * rho) == 1 / n;
        end
        
        
        % Entangled States
        
        function rep = getPauliRepresentation(Quantum, rho)
            Pauli = Quantum.getPauliMatrices();
            C = Pauli()
        end
        
        function res = isEntangled2(Quantum, psi)
            res = psi(1) * psi(4) ~= psi(2) * psi(3);
        end
    end
end