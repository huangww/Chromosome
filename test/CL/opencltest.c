#include <stdio.h>
#include <OpenCL/opencl.h>
#include "test.cl.h"
#define NUM_VALUES 1024

static void print_device_info(cl_device_id device)
{
	char name[128];
	char vendor[128];

	clGetDeviceInfo(device, CL_DEVICE_NAME, 128, name, NULL);
	clGetDeviceInfo(device, CL_DEVICE_VENDOR, 128, vendor, NULL);

}

static int validate(cl_float* input, cl_float* output)
{
	for (int i = 0; i < NUM_VALUES; ++i)
	{
		if (output[i] != (input[i] * input[i]))
		{
			fprintf(stdout, "Error: Element %d did not match experted output.\n", i);
			fprintf(stdout, 
					"	Saw %1.4f, expected %1.4f\n", output[i], input[i]*input[i]);
			fflush(stdout);
			return 0;
		}
	}
	return 1;
}

int main(int argc, char *argv[])
{
	char name[128];

	dispatch_queue_t queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_GPU, NULL);

	if (queue == NULL)
	{
		queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_CPU, NULL);
	}
	
	cl_device_id gpu = gcl_get_device_id_with_dispatch_queue(queue);
	clGetDeviceInfo(gpu, CL_DEVICE_NAME, 128, name, NULL);
	fprintf(stdout, "Created a dispatch queue using the %s\n", name);


	float* test_in = (float*)malloc(sizeof(cl_float) * NUM_VALUES);
	for (int i = 0; i < NUM_VALUES; ++i)
	{
		test_in[i] = (cl_float)i;
	}
	float* test_out = (float*)malloc(sizeof(cl_float) * NUM_VALUES);


	void* mem_in = gcl_malloc(sizeof(cl_float) * NUM_VALUES, test_in, 
			CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
	void* mem_out = gcl_malloc(sizeof(cl_float) * NUM_VALUES, NULL, 
			CL_MEM_WRITE_ONLY);

	dispatch_sync(queue, ^{
			size_t wgs;
			gcl_get_kernel_block_workgroup_info(square_kernel,
				CL_KERNEL_WORK_GROUP_SIZE,
				sizeof(wgs), &wgs, NULL);
			cl_ndrange range = {1, {0,0,0}, {NUM_VALUES, 0, 0}, {wgs, 0, 0}};
			square_kernel(&range, (cl_float*)mem_in, (cl_float*)mem_out);
			gcl_memcpy(test_out, mem_out, sizeof(cl_float) * NUM_VALUES);
			}); 
	if (validate(test_in, test_out))
	{
		fprintf(stdout, "All values were properly squared. \n");
	}

	gcl_free(mem_in);
	gcl_free(mem_out);

	free(test_in);
	free(test_out);

	dispatch_release(queue);

	return 0;
}

