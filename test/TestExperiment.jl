ret = Experiment(interval(35)/10, interval(1.0), interval(0.1), 128)

@test abs(ret.lo-0.196215) < (10)^(-3)