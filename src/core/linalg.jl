#=
Quaternion, DCM, and related linear-algebra utilities for attitude and
3-D geometry operations.
=#

"""
    sat_norm(n) -> Number

Saturate a near-unit scalar into ``[-1, 1]`` to guard against floating-point
error before inverse trigonometric calls.

# Arguments
- `n`: scalar cosine-like value that may lie slightly outside ``[-1, 1]``

# Returns
- saturated value in ``[-1, 1]`` when `|n|` is within `1 + 1e-6` of unity; otherwise `n` unchanged
"""
sat_norm(n) = abs(n) < 1 + 1e-6 ? min(n,max(n,-1),1) : n

"""
    shortest_arc_rotation(u, v) -> Vector

Compute the unit quaternion for the shortest-arc rotation taking vector `u`
onto vector `v`.

# Arguments
- `u`: source 3-vector
- `v`: destination 3-vector

# Returns
- scalar-first unit quaternion ``[q0, q1, q2, q3]``; identity `[1,0,0,0]` if `u ≈ v` for numerical stability
"""
function shortest_arc_rotation(u,v)
    θ = acos(sat_norm(dot(u,v)/(norm(u)*norm(v))))
    n = cross(u,v)
    n /= norm(n)
    if u ≈ v
        return [1,0,0,0]
    else
        return aa_to_quat(θ,n)
    end
end

"""
    aa_to_quat(angle, axis) -> Vector

Convert an axis-angle pair to a scalar-first unit quaternion.

# Arguments
- `angle`: rotation angle ``[rad]``
- `axis`: unit rotation axis (3-vector)

# Returns
- scalar-first unit quaternion ``[q0, q1, q2, q3]``
"""
function aa_to_quat(angle, axis)
    return [cos(angle/2), sin(angle/2)*axis...]
end

"""
    quat_to_aa(quat) -> (θ, n)

Convert a scalar-first unit quaternion to axis-angle form.

# Arguments
- `quat`: quaternion ``[q0, q1, q2, q3]``

# Returns
- `θ`: rotation angle ``[rad]``
- `n`: unit rotation axis (arbitrary `[1,0,0]` if `θ ≈ 0`)
"""
function quat_to_aa(quat)
    θ = 2*acos(quat[1])
    if θ ≈ 0
        n = [1,0,0] # axis undefined, choose arbitrarily
    else
        n = quat[2:4] / sin(θ/2)
    end
    return θ,n
end

"""
    quat_to_dcm(quat) -> Matrix

Convert a scalar-first unit quaternion to a direction cosine matrix (DCM).

# Arguments
- `quat`: quaternion ``[q0, q1, q2, q3]``

# Returns
- ``3 \\times 3`` rotation matrix mapping body to inertial (or vice versa per convention used here)
"""
function quat_to_dcm(quat)
    q0,q1,q2,q3 = quat
    return [
        1-2*(q2^2+q3^2) 2*(q1*q2+q0*q3) 2*(q1*q3-q0*q2);
        2*(q1*q2-q0*q3) 1-2*(q1^2+q3^2) 2*(q2*q3+q0*q1);
        2*(q1*q3+q0*q2) 2*(q2*q3-q0*q1) 1-2*(q1^2+q2^2)
    ]
end

"""
    dcm_to_quat(dcm) -> Vector

Convert a DCM to a scalar-first unit quaternion using a numerically stable
method (MATLAB-style branching on the matrix trace/diagonal).

# Arguments
- `dcm`: ``3 \\times 3`` direction cosine matrix

# Returns
- quaternion ``[q0, q1, q2, q3]``
"""
function dcm_to_quat(dcm)
    # See: https://www.mathworks.com/matlabcentral/answers/164746-bug-in-dcm2quat-function
    # Numerically-stable method as implemented by MATLAB!
    C11,C12,C13 = dcm[1,:]
    C21,C22,C23 = dcm[2,:]
    C31,C32,C33 = dcm[3,:]

    trace = tr(dcm)
    d = diag(dcm)

    if (trace > 0)
        sqtrp1 = sqrt( trace + 1.0 )
        
        q0 = 0.5*sqtrp1 
        q1 = (C23 - C32)/(2.0*sqtrp1)
        q2 = (C31 - C13)/(2.0*sqtrp1) 
        q3 = (C12 - C21)/(2.0*sqtrp1) 
    else
        if ((d[2] >= d[1]) && (d[2] >= d[3]))
            if C31 - C13 >= 0
                sqdip1 =  sqrt(d[2] - d[1] - d[3] + 1.0 )
            else
                sqdip1 = -sqrt(d[2] - d[1] - d[3] + 1.0 )
            end
            
            q2 = 0.5*sqdip1 
            
            if ( sqdip1 != 0 )
                sqdip1 = 0.5/sqdip1
            end
            
            q0 = (C31 - C13)*sqdip1 
            q1 = (C12 + C21)*sqdip1 
            q3 = (C23 + C32)*sqdip1 
        elseif (d[3] >= d[1])
            if C12 - C21 >= 0
                sqdip1 =  sqrt(d[3] - d[1] - d[2] + 1.0 )
            else
                sqdip1 = -sqrt(d[3] - d[1] - d[2] + 1.0 )
            end
            
            q3 = 0.5*sqdip1 
            
            if ( sqdip1 != 0 )
                sqdip1 = 0.5/sqdip1
            end
            
            q0 = (C12 - C21)*sqdip1
            q1 = (C31 + C13)*sqdip1 
            q2 = (C23 + C32)*sqdip1 
        else
            if C12 - C21 >= 0
                sqdip1 =  sqrt(d[1] - d[2] - d[3] + 1.0 )
            else
                sqdip1 = -sqrt(d[1] - d[2] - d[3] + 1.0 )
            end

            q1 = 0.5*sqdip1 
            
            if ( sqdip1 != 0 )
                sqdip1 = 0.5/sqdip1
            end
            
            q0 = (C23 - C32)*sqdip1 
            q2 = (C12 + C21)*sqdip1 
            q3 = (C31 + C13)*sqdip1 
        end
    end

    return [q0,q1,q2,q3]
end

"""
    quat_slerp(q1, q2, f) -> Vector

Spherical linear interpolation between unit quaternions.

# Arguments
- `q1`: start quaternion
- `q2`: end quaternion
- `f`: interpolation fraction in ``[0, 1]`` (`0` → `q1`, `1` → `q2`)

# Returns
- interpolated unit quaternion at fraction `f` (equals `q1` if the angle is ≈ 0)
"""
function quat_slerp(q1,q2,f)
    Ω = acos(dot(q1,q2))
    if Ω ≈ 0
        slerp = q1
    else
        slerp = sin((1-f)*Ω) / sin(Ω) * q1 + sin(f*Ω) / sin(Ω) * q2
    end
    return slerp
end

"""
    skew(vec) -> Matrix

Form the ``3 \\times 3`` skew-symmetric matrix associated with a 3-vector.

# Arguments
- `vec`: 3-vector ``[x, y, z]``

# Returns
- skew-symmetric matrix ``[v]_×`` such that ``[v]_× w = v \\times w``
"""
function skew(vec)
    x,y,z = vec
    return [
         0 -z  y;
         z  0 -x;
        -y  x  0
    ]
end

"""
    skew_projR4(vec) -> Matrix

Form the ``4 \\times 4`` quaternion kinematic skew matrix associated with a
3-vector angular rate.

# Arguments
- `vec`: 3-vector angular-rate components ``[x, y, z]``

# Returns
- ``4 \\times 4`` matrix used in quaternion kinematics ``\\dot{q} = \\tfrac{1}{2} Ω(ω) q``
"""
function skew_projR4(vec)
    x,y,z = vec
    return [
         0 -x -y -z;
         x  0  z -y;
         y -z  0  x;
         z  y -x  0
    ]
end
