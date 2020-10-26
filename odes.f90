      LOSS = -RATE(4)*Y(6)*D-RATE(5)*Y(3)*D-RATE(9)-RATE(12)-RATE(14)&
    &*Y(6)*D-RATE(17)*Y(4)*D
      PROD = +RATE(1)*Y(6)*Y(3)*D+2*RATE(2)*Y(3)*Y(3)*D+2&
    &*RATE(3)*Y(3)+2*RATE(4)*Y(1)*Y(6)*D+3*RATE(5)*Y(1)*Y(3)*D+RATE(7)&
    &*Y(3)+2*RATE(8)*Y(3)+RATE(11)*Y(6)+RATE(13)*Y(4)*Y(3)*D+RATE(16)&
    &*Y(6)+RATE(19)*Y(2)
      YDOT(1) = PROD+Y(1)*LOSS
      LOSS = -RATE(19)
      PROD = +RATE(7)*Y(3)+RATE(9)*Y(1)+RATE(12)*Y(1)
      YDOT(2) = PROD+Y(2)*LOSS
      LOSS = -RATE(1)*Y(6)*D-2*RATE(2)*Y(3)*D-RATE(3)-RATE(5)&
    &*Y(1)*D-RATE(7)-RATE(8)-RATE(13)*Y(4)*D
      PROD = +RATE(1)*Y(6)*Y(3)*D+RATE(2)*Y(3)*Y(3)*D+RATE(14)&
    &*Y(1)*Y(6)*D
      YDOT(3) = PROD+Y(3)*LOSS
      LOSS = -RATE(6)-RATE(10)-RATE(13)*Y(3)*D-RATE(15)-RATE(17)*Y(1)*D
      PROD = +RATE(1)*Y(6)*Y(3)*D+RATE(4)*Y(1)*Y(6)*D+RATE(11)&
    &*Y(6)+RATE(14)*Y(1)*Y(6)*D+RATE(16)*Y(6)+RATE(18)*Y(5)
      YDOT(4) = PROD+Y(4)*LOSS
      LOSS = -RATE(18)
      PROD = +RATE(6)*Y(4)+RATE(10)*Y(4)+RATE(15)*Y(4)
      YDOT(5) = PROD+Y(5)*LOSS
      
      LOSS = -RATE(1)*Y(3)*D-RATE(4)*Y(1)*D-RATE(11)-RATE(14)*Y(1)*D-RATE(16) ! what is D

      PROD = +RATE(13)*Y(4)*Y(3)*D+RATE(17)*Y(4)*Y(1)*D
      
      YDOT(6) = PROD+Y(6)*LOSS
