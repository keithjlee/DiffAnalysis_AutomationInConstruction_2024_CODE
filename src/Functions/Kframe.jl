"""
    k_frame(E::Float64, G::Float64, A::Float64, L::Float64, Ix::Float64, Iy::Float64, J::Float64)

Element stiffness matrix for a fixed-fixed frame element in LCS
"""
function k_frame(E::Float64, G::Float64, A::Float64, L::Float64, Ix::Float64, Iy::Float64, J::Float64)

    E / L^3 * [
        A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 12Ix 0 0 0 6L*Ix 0 -12Ix 0 0 0 6L*Ix;
        0 0 12Iy 0 -6L*Iy 0 0 0 -12Iy 0 -6L*Iy 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 -6L*Iy 0 4L^2*Iy 0 0 0 6L*Iy 0 2L^2*Iy 0;
        0 6L*Ix 0 0 0 4L^2*Ix 0 -6L*Ix 0 0 0 2L^2*Ix;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -12Ix 0 0 0 -6L*Ix 0 12Ix 0 0 0 -6L*Ix;
        0 0 -12Iy 0 6L*Iy 0 0 0 12Iy 0 6L*Iy 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 -6L*Iy 0 2L^2*Iy 0 0 0 6L*Iy 0 4L^2*Iy 0;
        0 6L*Ix 0 0 0 2L^2*Ix 0 -6L*Ix 0 0 0 4L^2*Ix
        ]
end

"""
    k_frame2(E::Float64, G::Float64, A::Float64, L::Float64, Ix::Float64, Iy::Float64, J::Float64)

Element stiffness matrix for a fixed-fixed frame element in LCS. experimental version for explicit rrule testing
"""
function _k_frame(E::Float64, G::Float64, A::Float64, L::Float64, Ix::Float64, Iy::Float64, J::Float64)

    return [
        E*A/L 0 0 0 0 0 -A*E/L 0 0 0 0 0;
        0 12Ix*E/L^3 0 0 0 6Ix*E/L^2 0 -12Ix*E/L^3 0 0 0 6Ix*E/L^2;
        0 0 12Iy*E/L^3 0 -6Iy*E/L^2 0 0 0 -12Iy*E/L^3 0 -6Iy*E/L^2 0;
        0 0 0 G*J/L 0 0 0 0 0 -G*J/L 0 0;
        0 0 -6Iy*E/L^2 0 4Iy*E/L 0 0 0 6Iy*E/L^2 0 2Iy*E/L 0;
        0 6Ix*E/L^2 0 0 0 4Ix*E/L 0 -6Ix*E/L^2 0 0 0 2Ix*E/L;
        -A*E/L 0 0 0 0 0 A*E/L 0 0 0 0 0;
        0 -12Ix*E/L^3 0 0 0 -6Ix*E/L^2 0 12Ix*E/L^3 0 0 0 -6Ix*E/L^2;
        0 0 -12Iy*E/L^3 0 6Iy*E/L^2 0 0 0 12Iy*E/L^3 0 6Iy*E/L^2 0;
        0 0 0 -G*J/L 0 0 0 0 0 G*J/L 0 0;
        0 0 -6Iy*E/L^2 0 2Iy*E/L 0 0 0 6Iy*E/L^2 0 4Iy*E/L 0;
        0 6Ix*E/L^2 0 0 0 2Ix*E/L 0 -6Ix*E/L^2 0 0 0 4Ix*E/L
        ]

end

function ChainRulesCore.rrule(::typeof(_k_frame), E::Float64, G::Float64, A::Float64, L::Float64, Ix::Float64, Iy::Float64, J::Float64)

    k = _k_frame(E, G, A, L, Ix, Iy, J)

    function _k_frame_pullback(k̄)

        dA = [
            E/L 0 0 0 0 0 -E/L 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            -E/L 0 0 0 0 0 E/L 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            ]

        ∇A = dot(k̄, dA)

        dL = [
            -E*A/L^2 0 0 0 0 0 A*E/L^2 0 0 0 0 0;
            0 -36Ix*E/L^4 0 0 0 -12Ix*E/L^3 0 36Ix*E/L^4 0 0 0 -12Ix*E/L^3;
            0 0 -36Iy*E/L^4 0 12Iy*E/L^3 0 0 0 36Iy*E/L^4 0 12Iy*E/L^3 0;
            0 0 0 -G*J/L^2 0 0 0 0 0 G*J/L^2 0 0;
            0 0 12Iy*E/L^3 0 -4Iy*E/L^2 0 0 0 -12Iy*E/L^3 0 -2Iy*E/L^2 0;
            0 -12Ix*E/L^3 0 0 0 -4Ix*E/L^2 0 12Ix*E/L^3 0 0 0 -2Ix*E/L^2;
            A*E/L^2 0 0 0 0 0 -A*E/L^2 0 0 0 0 0;
            0 36Ix*E/L^4 0 0 0 12Ix*E/L^3 0 -36Ix*E/L^4 0 0 0 12Ix*E/L^3;
            0 0 36Iy*E/L^4 0 -12Iy*E/L^3 0 0 0 -36Iy*E/L^4 0 -12Iy*E/L^3 0;
            0 0 0 G*J/L^2 0 0 0 0 0 -G*J/L^2 0 0;
            0 0 12Iy*E/L^3 0 -2Iy*E/L^2 0 0 0 -12Iy*E/L^3 0 -4Iy*E/L^2 0;
            0 -12Ix*E/L^3 0 0 0 -2Ix*E/L^2 0 12Ix*E/L^3 0 0 0 -4Ix*E/L^2
            ]

        ∇L = dot(k̄, dL)

        dIx = [
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 12E/L^3 0 0 0 6E/L^2 0 -12E/L^3 0 0 0 6E/L^2;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 6E/L^2 0 0 0 4E/L 0 -6E/L^2 0 0 0 2E/L;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 -12E/L^3 0 0 0 -6E/L^2 0 12E/L^3 0 0 0 -6E/L^2;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 6E/L^2 0 0 0 2E/L 0 -6E/L^2 0 0 0 4E/L
            ]

        ∇Ix = dot(k̄, dIx)

        dIy = [
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 12E/L^3 0 -6E/L^2 0 0 0 -12E/L^3 0 -6E/L^2 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 -6E/L^2 0 4E/L 0 0 0 6E/L^2 0 2E/L 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 -12E/L^3 0 6E/L^2 0 0 0 12E/L^3 0 6E/L^2 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 -6E/L^2 0 2E/L 0 0 0 6E/L^2 0 4E/L 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            ]

        ∇Iy = dot(k̄, dIy)

        dJ = [
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 G/L 0 0 0 0 0 -G/L 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 -G/L 0 0 0 0 0 G/L 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;
            ]

        ∇J = dot(k̄, dJ)

        return (
            NoTangent(),
            NoTangent(),
            NoTangent(),
            ∇A,
            ∇L,
            ∇Ix,
            ∇Iy,
            ∇J
        )
    end

    return k, _k_frame_pullback

end

function k(::Type{Asap.FixedFixed}, E, G, A, L, Ix, Iy, J)

    return E / L^3 * [
        A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 12Ix 0 0 0 6L*Ix 0 -12Ix 0 0 0 6L*Ix;
        0 0 12Iy 0 -6L*Iy 0 0 0 -12Iy 0 -6L*Iy 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 -6L*Iy 0 4L^2*Iy 0 0 0 6L*Iy 0 2L^2*Iy 0;
        0 6L*Ix 0 0 0 4L^2*Ix 0 -6L*Ix 0 0 0 2L^2*Ix;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -12Ix 0 0 0 -6L*Ix 0 12Ix 0 0 0 -6L*Ix;
        0 0 -12Iy 0 6L*Iy 0 0 0 12Iy 0 6L*Iy 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 -6L*Iy 0 2L^2*Iy 0 0 0 6L*Iy 0 4L^2*Iy 0;
        0 6L*Ix 0 0 0 2L^2*Ix 0 -6L*Ix 0 0 0 4L^2*Ix
        ]

end

function k(::Type{Asap.FreeFixed}, E, G, A, L, Ix, Iy, J)

    return E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Ix 0 0 0 0 0 -3Ix 0 0 0 3L*Ix;
        0 0 3Iy 0 0 0 0 0 -3Iy 0 -3L*Iy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Ix 0 0 0 0 0 3Ix 0 0 0 -3L*Ix;
        0 0 -3Iy 0 0 0 0 0 3Iy 0 3L*Iy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iy 0 0 0 0 0 3L*Iy 0 3L^2*Iy 0;
        0 3L*Ix 0 0 0 0 0 -3L*Ix 0 0 0 3L^2*Ix    
        ]

end

function k(::Type{Asap.FixedFree}, E, G, A, L, Ix, Iy, J)

    return E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Ix 0 0 0 3L*Ix 0 -3Ix 0 0 0 0;
        0 0 3Iy 0 -3L*Iy 0 0 0 -3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iy 0 3L^2*Iy 0 0 0 3L*Iy 0 0 0;
        0 3L*Ix 0 0 0 3L^2*Ix 0 -3L*Ix 0 0 0 0 ;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Ix 0 0 0 -3L*Ix 0 3Ix 0 0 0 0;
        0 0 -3Iy 0 3L*Iy 0 0 0 3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
        ]

end

function k(::Type{Asap.Joist}, E, G, A, L, Ix, Iy, J)

    return E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Ix 0 0 0 3L*Ix 0 -3Ix 0 0 0 0;
        0 0 3Iy 0 -3L*Iy 0 0 0 -3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iy 0 3L^2*Iy 0 0 0 3L*Iy 0 0 0;
        0 3L*Ix 0 0 0 3L^2*Ix 0 -3L*Ix 0 0 0 0 ;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Ix 0 0 0 -3L*Ix 0 3Ix 0 0 0 0;
        0 0 -3Iy 0 3L*Iy 0 0 0 3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
        ]

end