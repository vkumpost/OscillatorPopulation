
err = OscillatorPopulationError("test message!!!")
@test err isa Exception
@test sprint(showerror, err) == "OscillatorPopulationError: test message!!!"
