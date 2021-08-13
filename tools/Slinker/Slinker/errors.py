#=======================================================================================================================
#
#   Slinker - Error handling
#   Author: Breon Schmidt
#   License: MIT
#
#=======================================================================================================================

''' --------------------------------------------------------------------------------------------------------------------
Imports
---------------------------------------------------------------------------------------------------------------------'''

import sys

''' --------------------------------------------------------------------------------------------------------------------
Functions
---------------------------------------------------------------------------------------------------------------------'''

def error(code, support=False):

	'''
		Slinker errors are encoded in a three digit code:
		first digit - 0: Warning, 1: Error (terminating)
		second - 0/1 - 0: Missing, 1: Incorrect
		third - 0/1 - 0: Argument, 1: Key, 2: Type
	'''

	error = str(code)
	message = [error+" -"]
	exit = False

	'''First Digit'''
	if error[0] == "0":
		message.append("WARNING:")
	elif error[0] == "1":
		exit = True
		message.append("ERROR:")

	'''Second Digit'''
	message.append(support)
	if error[1] == "0":
		message.append("is a missing")
	if error[1] == "1":
		message.append("is an incorrect")

	'''Third Digit'''
	if error[2] == "0":
		message.append("argument.")
	if error[2] == "1":
		message.append("key.")
	if error[2] == "2":
		message.append("type.")

	'''Send message to user of the bad news'''
	if exit:
		message.append("TERMINATING.")

	print(" ".join(message))

	if exit:
		sys.exit()





