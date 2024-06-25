sat_norm(n) = abs(n) < 1 + 1e-6 ? min(n,max(n,-1),1) : n # saturate normalized input between -1 and 1 if slightly imprecise due to floating point error

function shortest_arc_rotation(u,v)
    # Shortest-arc rotation between two vectors u,v
    θ = acos(sat_norm(dot(u,v)/(norm(u)*norm(v))))
    n = cross(u,v)
    n /= norm(n)
    if u ≈ v
        return [1,0,0,0]
    else
        return aa_to_quat(θ,n)
    end
end

function aa_to_quat(angle, axis)
    return [cos(angle/2), sin(angle/2)*axis...]
end

function quat_to_aa(quat)
    θ = 2*acos(quat[1])
    if θ ≈ 0
        n = [1,0,0] # axis undefined, choose arbitrarily
    else
        n = quat[2:4] / sin(θ/2)
    end
    return θ,n
end

function quat_to_dcm(quat)
    q0,q1,q2,q3 = quat
    return [
        1-2*(q2^2+q3^2) 2*(q1*q2+q0*q3) 2*(q1*q3-q0*q2);
        2*(q1*q2-q0*q3) 1-2*(q1^2+q3^2) 2*(q2*q3+q0*q1);
        2*(q1*q3+q0*q2) 2*(q2*q3-q0*q1) 1-2*(q1^2+q2^2)
    ]
end

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

function quat_slerp(q1,q2,f)
    Ω = acos(dot(q1,q2))
    if Ω ≈ 0
        slerp = q1
    else
        slerp = sin((1-f)*Ω) / sin(Ω) * q1 + sin(f*Ω) / sin(Ω) * q2
    end
    return slerp
end

function skew(vec)
    x,y,z = vec
    return [
         0 -z  y;
         z  0 -x;
        -y  x  0
    ]
end

function skew_projR4(vec)
    x,y,z = vec
    return [
         0 -x -y -z;
         x  0  z -y;
         y -z  0  x;
         z  y -x  0
    ]
end