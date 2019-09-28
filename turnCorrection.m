function err = turnCorrection(theta,theta_d)
    if (theta>0 && theta_d<0 && abs(theta-theta_d)>pi)
        err = 2*pi-abs(theta-theta_d);
    elseif (theta<0 && theta_d>0 && abs(theta-theta_d)>pi)
        err = -(2*pi-abs(theta-theta_d));
    else
        err = theta-theta_d;
    end
end
