#inlucde planet.hpp

Planet::planet(double m, vec position){
  mass = m;
  pos = position;
  vel = velocity;
}



Planet::distanceOther(){
  vec<double> dr = pos - otherPlanet.pos;
  return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));
}

Planet::gravitationalForce(Planet otherPlanet){
  double r = distanceOther(Planet otherPlanet);
  vec<double> Fg = G * mass * otherPlanet.mass / (r * r)
  return Fg
}
