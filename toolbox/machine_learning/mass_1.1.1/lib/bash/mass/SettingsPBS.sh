#!/bin/bash

setSubmitCMD()
{
	# Getting the arguments
	_memory=$1
	shift
	_holdids="$*"

	# Getting the options for PBS scheduler into variables
	# This is done so that the user can adapt these options to match them with his own scheduler
	_OptExec="qsub"
	_OptPriority="-p"
	_OptMemory="-l mem=${_memory}gb,walltime=4:00:00"
	_OptTerse=""
	_OptJoin="-j oe"
	_OptOutput="-o ${log}\$PBS_JOBNAME-\$PBS_JOBID.log"
	_OptHoldID="-W depend=afterany:"


	# Setting up the basic command and then the other options will be added onto it
	submitCMD="${_OptExec} ${_OptTerse} ${_OptJoin} ${_OptOutput} ${_OptMemory}"


	# check if priority variable
	if [ -n "${priority}" ]
	then
		submitCMD="${submitCMD} ${_OptPriority} ${priority}";
	fi

	# check if the user provided a hold jobid
	# Modify this part of the code depending on how the job ids need to be appended
	if [ -n "${_holdids}" ]
	then
		_holdidsString=`echo ${_holdids} | sed 's/ /:/g'`
		submitCMD="${submitCMD} ${_OptHoldID}${_holdidsString}"
	fi

	echo $submitCMD
}

setDeleteJob()
{
	# Getting the job IDs
	_jobid="$*"

	# Getting the options for PBS scheduler into variables
	# This is done so that the user can adapt these options to match them with his own scheduler
	_OptJobDeleteExec="qdel"					# Executable name for deleting jobs

	echo ${_OptJobDeleteExec}
}
