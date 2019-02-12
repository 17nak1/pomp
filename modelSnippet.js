
snippet = {}
let mathLib = require('./mathLib')

snippet.skeleton = function (params, t, N, pop, birthrate) {
  var seas, dy = []
  var R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] 
  var beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma
  var S = N[0], E = N[1], R = N[2], I = N[3]
  var va = snippet.skeletonVaccine(t)
  var tt = (t - Math.floor(t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - amplitude
  }
  var Beta = beta0 * seas / pop
  dy[0] = birthrate * (1 - va) - Beta * S * I - mu * S
  dy[1] = Beta * S * I - (sigma + mu) * E
  dy[2] = gamma * I - mu * R + birthrate * va
  dy[3] = sigma * E - (gamma + mu) * I
  dy[4] = gamma * I
  return dy
}

snippet.skeletonVaccine = function(t) {
  var vaccineRate
  if (t < 1968)
    vaccineRate = 0
  else if (t >= 1968 && t <= 1969)
    vaccineRate = 0.33
  else if (t >= 1969 && t <= 1970)
    vaccineRate = 0.46
  else if (t >= 1970 && t <= 1971)
    vaccineRate = 0.51
  else if (t >= 1971 && t <= 1972)
    vaccineRate = 0.53
  else if (t >= 1972 && t <= 1973)
    vaccineRate = 0.52
  else if (t >= 1973 && t <= 1974)
    vaccineRate = 0.46
  else if (t >= 1974 && t <= 1975)
    vaccineRate = 0.46
  else if (t >= 1975 && t <= 1976)
    vaccineRate = 0.48
  else if (t >= 1976 && t <= 1977)
    vaccineRate = 0.48
  else if (t >= 1977 && t <= 1978)
    vaccineRate = 0.51
  else if (t >= 1978 && t <= 1979)
    vaccineRate = 0.53;
  else if (t >= 1979 && t <= 1980)
    vaccineRate = 0.55;
  else if (t >= 1980 && t <= 1981)
    vaccineRate = 0.58;
  else if (t >= 1981 && t <= 1982)
    vaccineRate = 0.60
  else if (t >= 1982 && t <= 1983)
    vaccineRate = 0.63
  else if (t >= 1983 && t <= 1984)
    vaccineRate = 0.68
  else if (t >= 1984 && t <= 1985)
    vaccineRate = 0.71
  else if (t >= 1985 && t <= 1988)
    vaccineRate = 0.76
  else if (t >= 1988 && t <= 1989)
    vaccineRate = 0.814
  else if (t >= 1989 && t <= 1990)
    vaccineRate = 0.9488
  else if (t >= 1990 && t <= 1991)
    vaccineRate = 0.9818
  else if (t >= 1991 && t <= 1992)
    vaccineRate = 0.90
  else if (t >= 1992 && t <= 1993)
    vaccineRate = 0.92
  else if (t >= 1993 && t <= 1994)
    vaccineRate = 0.91
  else if (t >= 1994 && t <= 1995)
    vaccineRate = 0.91
  else if (t >= 1995 && t <= 1996)
    vaccineRate = 0.92
  else if (t >= 1996 && t <= 1997)
    vaccineRate = 0.92
  else if (t >= 1997 && t <= 1998)
    vaccineRate = 0.91
  else if (t >= 1998 && t <= 1999)
    vaccineRate = 0.88
  else if (t >= 1999 && t <= 2000)
    vaccineRate = 0.88
  else if (t >= 2000 && t <= 2001)
    vaccineRate = 0.87
  else if (t >= 2001 && t <= 2002)
    vaccineRate = 0.84
  else if (t >= 2002 && t <= 2003)
    vaccineRate = 0.82
  else if (t >= 2003 && t <= 2004)
    vaccineRate = 0.80
  else if (t >= 2004 && t <= 2005)
    vaccineRate = 0.81
  else if (t >= 2005 && t <= 2006)
    vaccineRate = 0.84
  else if (t >= 2006 && t <= 2007)
    vaccineRate = 0.85
  else if (t >= 2007 && t <= 2008)
    vaccineRate = 0.85
  else if (t >= 2008 && t <= 2009)
    vaccineRate = 0.85
  else if (t >= 2009 && t <= 2010)
    vaccineRate = 0.88
  else
    vaccineRate = 0.89
  return vaccineRate
}
snippet.initz = function(pop, S, E, R, I) {
  var m = pop / (S + E + R + I),
    S = Math.round(m * S),
    E = Math.round(m * E),
    R = Math.round(m * R),
    I = Math.round(m * I),
    H = 0
  return [S, E, R, I, H]
}

snippet.dmeasure = function (rho, psi, H, dCases, giveLog = 1) {
  var lik
  var mn = rho * H
  var v = mn * (1.0 - rho + psi * psi * mn)
  var tol = 1.0e-18
  var modelCases = Number(dCases)
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
    }
    if (giveLog) lik = Math.log(lik)
  } else {
    lik = (giveLog) ? 0 : 1;
  }
  return lik
}

snippet.rmeasure = function (H, rho, psi) {
  var mn = rho * H
  var v = mn * (1.0 - rho + psi * psi * mn)
  var tol = 1.0e-18
  var cases = mathLib.rnorm(mn, Math.sqrt(v) + tol)
  if (cases > 0) {
    cases = Math.round(cases)
  } else {
    cases = 0
  }
  return cases
}

module.exports = snippet