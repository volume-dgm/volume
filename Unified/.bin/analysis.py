processesCount = 32

globalStepTime = 0
stepsCount = 0

for p in range(processesCount): 
	with open('out/profiling' + '[' + str(p) + '].txt', 'r') as f:
		stepTimes = 0
		advancePhase = 0
		sync = 0
		allPhases = 0
		for line in f:
			if 'Step:' in line:
				stepTimes += float(line.split(':')[1])
				stepsCount += 1
			if 'Advance phase:' in line:
				advancePhase += float(line.split(':')[1])
			if 'Synchronization:' in line:
				sync += float(line.split(':')[1])

		globalStepTime += stepTimes

		print 'Process ', p, 'Quality ', float(advancePhase) / stepTimes, ' Sync ', sync / stepTimes 
		print '  Other ', (1 - (advancePhase + sync) / stepTimes) 

print ''
print 'Average step time: ', float(globalStepTime) / stepsCount