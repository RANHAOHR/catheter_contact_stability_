function g = se3rot(w,q,t)
% Map se(3) to SE(3) in case of pure rotation.
% From eq 2.40 MLS in [1].

    if size(w,2)>1 || size(q,2)>1
        disp('error: use row vector');
        return
    elseif size(w,1)~=3 || size(q,1)~=3
        disp('error: vector has length not equal to 3');
        return
    else
        % Find the rotation matrix
        R = so3rot(w,t);
        % Find g using eq 2.40
        g = [R, (eye(3)-R)*q; 0 0 0 1];
    end
    
end
