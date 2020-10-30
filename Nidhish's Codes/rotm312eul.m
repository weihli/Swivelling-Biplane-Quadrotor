function eul = rotm312eul(R) % Rotation matrix to 312 Euler angle, eul = (yaw,roll,pit)

eul(1,1) = atan2(-R(1,2),R(2,2)); % yaw
eul(2,1) = asin(R(3,2));  % roll
eul(3,1) = atan2(-R(3,1),R(3,3));  % pitch

end