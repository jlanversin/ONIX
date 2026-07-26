import onix.compute as compute

def pade(At, N):

	linalg = compute.get_expm_module()

	with compute.backend_device():
		At_backend = compute.asarray(At)
		N_backend = compute.asarray(N)
		expAt = linalg.expm(At_backend)
		result = expAt.dot(N_backend)
		compute.synchronize()

	return compute.to_numpy(result)
