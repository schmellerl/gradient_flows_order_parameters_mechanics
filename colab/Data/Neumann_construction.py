def Neumann_construction(theta_in_liquid, theta_in_solid, gamma_sa ):

  theta_in_liquid = theta_in_liquid*pi/180   
  theta_in_solid  = theta_in_liquid*pi/180

  theta_2 = 180*pi/180 - theta_outer 
  theta_1 = theta_inner - theta_2

  gamma_la = gamma_sa * 1/(cos(theta_1)+cos(theta_2)*sin(theta_1)/sin(theta_2) )
  gamma_sl = sin(theta_1)/sin(theta_2) * gamma_la

  gamma_s = (gamma_sa + gamma_sl - gamma_la)/2
  gamma_l = gamma_sl - gamma_s 
  gamma_a = gamma_sa - gamma_s
  
  return gamma_s, gamma_l, gamma_a

