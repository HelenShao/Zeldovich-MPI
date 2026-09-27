import numpy as np
import pyfftw
import os

# Steps for MPI (from zeldovich_mpi_driver.cpp) 
# * 7. MULTI-BATCH LOOP (for each batch) 
#  *    a. Get batch’s (y_primary, y_mirror)
#  *    b. Alloc local_y_slices + fft_stage_2d → generate + 2D FFT → free stage
#  *    c. calculate_batch_send_recv_counts; pack_slices_to_send_buffer; free local_y_slices
#  *    d. MPI_Alltoallv_c(send_buffer -> recv_buffer)
#  *    e. Update src_write_cursor; free per-batch send buffer

# Tell pyfftw to use the same number of threads as your system
pyfftw.interfaces.cache.enable()

#define dtype
record_dtype=np.dtype([
    ('x', np.int32),
    ('y', np.int32),
    ('z', np.int32),
    ('D_real', np.float64),
    ('D_imag', np.float64)
])

#parse parameter file for values
def parse_par2(path):
    """Parse a .par2 file into a dict of {key: value}.
    Handles 'key = value' lines, trailing/full-line '#' comments, and quoted
    strings. Values that look numeric are returned as float; otherwise as
    stripped strings.
    """
    params = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.split('#', 1)[0].strip()  # drop inline/full-line comments
            if not line or '=' not in line:
                continue
            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip().strip('"').strip("'")
            try:
                val = float(val)
            except ValueError:
                pass
            params[key] = val
    return params

#Obtain the PPD from data
sample_records = np.frombuffer(open("bin_files/Phi_gaussian_rank0_y0.bin", 'rb').read(), dtype=record_dtype)
N = int(round(np.sqrt(sample_records.shape[0]))) 

# Create aligned arrays for phi_gaussian(k) and |k|
phi_k = pyfftw.empty_aligned((N, N, N), dtype='complex128')
k_arr=pyfftw.empty_aligned((N,N,N), dtype='float64') 

#Get the BoxSize (L) from the parameter file in order to appropriately scale the |k| array
# BoxSize is read from the same param.par2 that run.sh feeds to Zeldovich_MPI,
# so this always matches the run that produced the Phi_gaussian_*.bin files.
PARAM_FILE = os.environ.get("PARAM_FILE", ".runtime/param.par2")
par2_params = parse_par2(PARAM_FILE)
L = par2_params["BoxSize"] ##think about units

#upload data from the dumped files into phi_k
for y in range(N//2 +1):
    filename = f"bin_files/Phi_gaussian_rank0_y{y}.bin"
    records = np.frombuffer(open(filename, 'rb').read(), dtype=record_dtype)
    phi_k[records['x'], records['y'], records['z']] = (records['D_real'] + 1j*records['D_imag']) 
    kfreq_1d = np.fft.fftfreq(N) * N          # [0, 1, ..., N/2-1, -N/2, ..., -1]
    KX, KY, KZ = np.meshgrid(kfreq_1d, kfreq_1d, kfreq_1d, indexing='ij')
    k_arr = ((2*np.pi) / L) * np.sqrt(KX**2 + KY**2 + KZ**2) 


##Handling Hermitian Symmetry##
index=np.arange(N) ##[ 0 1 2...N-1]
mirror= (-index) % N ##[0 N-1 N-2 ... 1]

y_index = np.arange(1,N//2) ##[1 2 ... N/2-1] excludes 250 if N=500
y_mirror = N-y_index ##[N-1 N-2 ... N/2 +1] excludes 250

phi_k[np.ix_(mirror,y_mirror,mirror)] = np.conj(phi_k[np.ix_(index, y_index, index)]) #This is Hermitian symmetry, it matches negative k values (indexed as -x --> N-x) to the conjugate of the positive k values

##Self-conjugate planes: y=0 and y=Nhalf need internal mirroring
for y0 in (0, N//2):
   z_index = np.arange(1, N//2)
   phi_k[np.ix_(mirror, [y0], N - z_index)] = np.conj(phi_k[np.ix_(index, [y0], z_index)]) #same mirroring process but excluding z=0 and z=N//2
   x_index = np.arange(1, N//2)
   phi_k[N - x_index, y0, 0]      = np.conj(phi_k[x_index, y0, 0]) #handling z=0 and z=N//2 cases (excluding x=0 and x=N//2)
   phi_k[N - x_index, y0, N//2]  = np.conj(phi_k[x_index, y0, N//2])
   for xx in (0, N//2):
       for zz in (0, N//2):
           phi_k[xx, y0, zz] = phi_k[xx, y0, zz].real  # force real at DC/Nyquist corners 

#In C, here we should free the index, mirror, y_index and y_mirror arrays



##Formulas for padding and cropping arrays (which we do to prevent aliasing)##
def pad_kspace(array, N, N2):
    pl = (N2 - N) // 2 # calculates how many zeros to pad on the left side of the array to get the array from N to N2
    pr = (N2 - N) - pl # calculates how many zeros to pad on the right side of the array
    out = np.fft.fftshift(array, axes=(0, 1, 2)) # shifts the array so that the back half of the data is now in the beginning entries (aka shifts the 0 mode to the center of the array) e.g. [0,1,2,3] --> [2,3,0,1]
    out = np.pad(out, ((pl, pr), (pl, pr), (pl, pr))) # pads the array with zeros on the left and right side of each axis
    return np.fft.ifftshift(out, axes=(0, 1, 2)) # shifts the array back so the padded zeros are now in the center

def crop_kspace(array_padded, N, N2):
    start = (N2 - N) // 2 # calculates how many zeros have been added to the left side get the original array from size N to size N2
    out = np.fft.fftshift(array_padded, axes=(0, 1, 2)) #shifts the array. Ensures the padded zeros are on the left and right ends of the array
    out = out[start:start+N, start:start+N, start:start+N] #only selects the center values
    return np.fft.ifftshift(out, axes=(0, 1, 2)) # shifts the array back to normal


##pad phi_k
phi_k_pad_data=pad_kspace(phi_k, N, 2*N) #created the padded phi array
#In C, here we should free phi_k (unpadded at this point) (not used again)
phi_k_pad=pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128')
phi_k_pad[:]=phi_k_pad_data #save padded phi_k into the pyfftw aligned empty array. Have to do it in these two steps to ensure the padded array is aligned for pyfftw.
#In C, free phi_k_pad_data (not used again)


##define function for convolution
#before using kernel function: use pyfftw empty for G, H and define G, H for your specific kernel
def kernel(G,H): #inputs G(k_1) and H(k_1) then returns your Q(k), which must then be multiplied by J(k) to get the final kernel 
    #pad the G and H arrays
    G_pad_data=pad_kspace(G, N, 2*N)
    #In C, free G
    G_pad=pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128')
    G_pad[:]=G_pad_data
    #In C, free G_pad_data (not used again)
    H_pad_data=pad_kspace(H, N, 2*N)
    #In C, free H
    H_pad=pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128')
    H_pad[:]=H_pad_data
    #In C, free H_pad_data (not used again)
    #Here, we're at 4 arrays in memory: phi_k_pad, k_array, G_pad, and H_pad

    Gphi_k = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128')
    Gphi_k=np.multiply(phi_k_pad,G_pad) #Multiply phi(k_1)*G(k_1)
    #In C, free G_pad


    Hphi_k = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128')
    Hphi_k=np.multiply(phi_k_pad, H_pad) #Multiply phi(k_2)*H(k_2)
    #In C, free H_pad
    
    #ffts of Gphi_k, and Hphi_k
    dummy = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128') #this takes us up to 5 arrays (phi_k_pad, k_arr, Gphi_k, Hphi_k, dummy)
    fft_xz = pyfftw.FFTW(dummy, dummy, 
                                axes=(2, 0), ##check with Helen's code which axis should go first
                                direction='FFTW_BACKWARD',
                                flags=('FFTW_MEASURE',))
    
    fft_y = pyfftw.FFTW(dummy, dummy, 
                        axes=(1,), 
                        direction='FFTW_BACKWARD',
                        flags=('FFTW_MEASURE',))

    #For multiple ranks, compute the xz fourier transformations in individual y-slices (individual y-slices are assigned to a given rank)
    fft_xz.update_arrays(Gphi_k, Gphi_k)   
    fft_xz()

    fft_xz.update_arrays(Hphi_k, Hphi_k)   
    fft_xz()
    #free dummy array
    #4 arrays: (phi_k_pad, k_arr, Gphi_k, Hphi_k)

    #here we will do AlltoAllv transformation to now divide up the xz segments into ranks 
    #I should also change it so here is where we enforce Hermitian symmetry on y-axis
    # This will take 2x the memory to transfer (currently at 4 total arrays, but only Gphi and Hphi would need to be transferred)
    #questions: what happends to phi_k_pad and k_arr during this transfer? Is it best to do this AlltoAllv transformation during every single kernel run? 
    # Or would it be best to format all the different Gphi and Hphi arrays for all the kernels (using much more memory) and then do the AlltoAllv transformation once?
    #If that is the best method, I would need to restructure the code a lot 
    
    #Now do the y fourier transformations xz segment by xz segment   
    Gphi_x = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128') #now we're up to 5 arrays
    fft_y.update_arrays(Gphi_k, Gphi_k)   
    fft_y()
    Gphi_x[:]= Gphi_k[:]   
    #In C, free Gphi_k 
    #back to 4 arrays

    Hphi_x = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128') #now we're up to 5 arrays
    fft_y.update_arrays(Hphi_k, Hphi_k)  
    fft_y()
    Hphi_x[:]= Hphi_k[:]
    #In C, free Hphi_k  
    #4 arrays: (phi_k_pad, k_arr, Gphi_x, Hphi_x) 


    #Q(x)=Gphi(x)*Hphi(x)
    Q_x_pad = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128') #5 arrays
    Q_x_pad=np.multiply(Gphi_x, Hphi_x)
    #In C, free Gphi_x and Hphi_x
    #3 arrays (phi_k_pad, k_arr, Q_x_pad)
    Q_x_pad *= (2*N)**3   # cancel the 1/Ntot_pad from pyfftw's normalized inverse transforms


    #fft of Q(x)
    dummy = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128') #this takes us up to 5 arrays (phi_k_pad, k_arr, Gphi_k, Hphi_k, dummy)
    fft_Q_y = pyfftw.FFTW(dummy, dummy, 
                axes=(1,), ##check with Helen's code which axis should go first
                direction='FFTW_FORWARD',
                flags=('FFTW_MEASURE',))

    #Here we will do another AlltoAllv to separate ranks by y-slice again
    
    fft_Q_xz = pyfftw.FFTW(dummy, dummy, 
                    axes=(2, 0), ##check with Helen's code which axis should go first
                    direction='FFTW_FORWARD',
                    flags=('FFTW_MEASURE',))

    
    #again, do the y fourier transformations xz segment by xz segment
    fft_Q_y.update_arrays(Q_x_pad, Q_x_pad)
    #clear dummy array from memory (not used again)
    fft_Q_y()

    #Here we do another AlltoAllv to separate ranks by y-slice again
    #do xz fourier transformation y-slice by y-slice
    Q_k_pad = pyfftw.empty_aligned((2*N, 2*N, 2*N), dtype='complex128') #4 arrays
    fft_Q_xz.update_arrays(Q_x_pad, Q_x_pad)
    fft_Q_xz()
    Q_k_pad[:]=Q_x_pad[:]
    #In C, free Q_x_pad
    #3 arrays: (phi_k_pad, k_arr, Q_k_pad)

    Q_k_pad[0][0][0]=0 ##Set DC mode as zero
    Q_k_data=crop_kspace(np.asarray(Q_k_pad),N,2*N)
    #In C, free Q_k_pad
    Q_k= pyfftw.empty_aligned((N, N, N), dtype='complex128')
    Q_k[:]=Q_k_data[:]
    #In C, free Q_k_data

    #3 arrays: phi_k_pad, k_arr, Q_k
    #phi_k_pad and k_arr are always around (must be used for all kernels), so after running a kernel function, we've only added one total array

    return Q_k


#Select which f_NL to calculate
local = par2_params["ZD_f_NL_local"] 
equil = par2_params["ZD_f_NL_equil"]
orthoCMB = par2_params["ZD_f_NL_orthoCMB"]


#f_NL local
if(local!=0):
    G_A1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    H_A1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    G_A1[:]=1
    H_A1[:]=1
    Q_A1=kernel(G_A1, H_A1)
    J_A1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    J_A1[:]=1
    K_A1=np.multiply(Q_A1,J_A1)
    #In C, free Q_A1, J_A1

    kernel_output=K_A1
    #In C, free K_A1


#input f_NL equil and f_NL orthoCMB (which use the same kernels)
if(equil!=0 or orthoCMB!=0):

    #Kernel K_A^I
    G_A1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    H_A1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    G_A1[:]=1
    H_A1[:]=1
    Q_A1=kernel(G_A1, H_A1)
    J_A1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    J_A1[:]=1
    K_A1=np.multiply(Q_A1,J_A1)
    #In C, free Q_A1, J_A1


    #Kernel K_C^I
    G_C1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    H_C1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    G_C1[:]=k_arr
    H_C1[:]=k_arr
    Q_C1=kernel(G_C1, H_C1)

    k_arr[0][0][0]=1 #momentarily modify the k_arr to avoid dividing by zero, which won't matter because we will set [0][0][0]=0 for the DC mode anyway
    J_C1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    J_C1=(1/3)*(1/np.multiply(k_arr,k_arr))
    J_C1[0][0][0]=0 #set DC mode to zero
    k_arr[0][0][0]=0 #set k_arr back to zero for DC mode
    K_C1=np.multiply(Q_C1,J_C1)
    #In C, free Q_C1, J_C1


    #Only one placement of the k weight is computed for K_B^I and K_B^II. Swapping G and H gives the same
    #real-space product (e.g. (d phi)*phi == phi*(d phi)), so the "second" kernels would equal the "first" ones,
    #and the coefficients below already refer to a single placement: 4*d^{-1}(phi d phi) + 2*nabla^{-2}(phi nabla^2 phi).
    #QuijotePNG paper's Table 4 writes K_B^I = (k_1+k_2)/k_3 (both placements, so i think using it with the coefficients 4 and 2
    #would unintetnionally double these terms?
    ##K_B^I (single placement: k_1/k_3)
    G_B1_first=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    H_B1_first=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    G_B1_first[:]=k_arr
    H_B1_first[:]=1
    Q_B1_first=kernel(G_B1_first,H_B1_first)

    J_B1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    k_arr[0][0][0]=1
    J_B1[:]=1/k_arr
    J_B1[0][0][0]=0
    k_arr[0][0][0]=0
    K_B1_first=np.multiply(Q_B1_first,J_B1)
    #In C, free Q_B1_first, J_B1

    # ##second kernel
    # G_B1_second=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    # H_B1_second=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    # G_B1_second[:]=1
    # H_B1_second[:]=k_arr
    # Q_B1_second=kernel(G_B1_second,H_B1_second)
    # J_B1=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    # k_arr[0][0][0]=1
    # J_B1[:]=1/k_arr
    # J_B1[0][0][0]=0
    # k_arr[0][0][0]=0
    # K_B1_second=np.multiply(Q_B1_second,J_B1)
    # #In C, free Q_B1_second, J_B1

    # K_B1=K_B1_first+K_B1_second

    #Kernel K_B^II is also a sum of two different kernels
    ##first kernel


    ##K_B^II (single placement: k_1^2/k_3^2)
    G_B2_first=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    H_B2_first=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    G_B2_first[:]=np.multiply(k_arr,k_arr)
    H_B2_first[:]=1
    Q_B2_first=kernel(G_B2_first,H_B2_first)

    J_B2=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    k_arr[0][0][0]=1
    J_B2[:]=1/np.multiply(k_arr, k_arr)
    J_B2[0][0][0]=0
    k_arr[0][0][0]=0
    K_B2_first=np.multiply(Q_B2_first,J_B2)
    #In C, free Q_B2_first, J_B2

    # ##second kernel
    # G_B2_second=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    # H_B2_second=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    # G_B2_second[:]=1
    # H_B2_second[:]=np.multiply(k_arr,k_arr)
    # Q_B2_second=kernel(G_B2_second,H_B2_second)

    # J_B2=pyfftw.empty_aligned((N,N,N), dtype='complex128')
    # k_arr[0][0][0]=1
    # J_B2[:]=1/np.multiply(k_arr, k_arr)
    # J_B2[0][0][0]=0
    # k_arr[0][0][0]=0
    # K_B2_second=np.multiply(Q_B2_second,J_B2)
    # #In C, free Q_B2_second, J_B2

    # K_B2=K_B2_first+K_B2_second

    #calculating phi_{NG}(k) for f_NL^{equil}
    if(equil!=0):
        # kernel_equil = -3*K_A1 - 6*K_C1 + 4*K_B1 + 2*K_B2
        kernel_equil = -3*K_A1 - 6*K_C1 + 4*K_B1_first + 2*K_B2_first
        kernel_output = kernel_equil
        #free K_A1, K_C1, K_B1_first, K_B2_first

    
    #calculating phi_{NG}(k) for f_NL^{orthoCMB}
    if(orthoCMB!=0):
        # kernel_orthoCMB = -9*K_A1 - 24*K_C1 + 10*K_B1 + 8*K_B2
        kernel_orthoCMB = -9*K_A1 - 24*K_C1 + 10*K_B1_first + 8*K_B2_first
        kernel_output = kernel_orthoCMB
        #free K_A1, K_C1, K_B1_first, K_B2_first


os.makedirs("bin_files", exist_ok=True)
os.makedirs("bin_files/bin_files_output", exist_ok=True)
 
x_all = np.arange(N)
z_all = np.arange(N)
# indexing='ij' with (z_all, x_all) gives arrays of shape (N, N) where axis 0
# is z and axis 1 is x; a C-order ravel() of these varies x fastest and z
# slowest, i.e. z=0 with all x-values, then z=1 with all x-values, etc.
zz, xx = np.meshgrid(z_all, x_all, indexing='ij')
#free x_all and z_all
xx_flat = xx.ravel()
zz_flat = zz.ravel()
#free xx and zz 

##upload final phi_non-gaussian(k) to files in y slices
for y in range(N//2 + 1):
    plane = kernel_output[:, y, :]              # shape (N, N), indexed [x, z]
    plane_zx = plane.T                          # reindex to [z, x] so ravel() matches xx_flat/zz_flat ordering
    kernel_real = plane_zx.real.ravel().astype(np.float64)
    kernel_imag = plane_zx.imag.ravel().astype(np.float64)
    yy_flat = np.full(xx_flat.shape, y, dtype=np.int32)
 
    out_records = np.zeros(xx_flat.shape[0], dtype=record_dtype)
    out_records['x'] = xx_flat.astype(np.int32)
    out_records['y'] = yy_flat
    out_records['z'] = zz_flat.astype(np.int32)
    out_records['D_real'] = kernel_real
    out_records['D_imag'] = kernel_imag
 
    filename = f"bin_files/bin_files_output/phi_output_rank0_y{y}.bin"
    with open(filename, 'wb') as f:
        f.write(out_records.tobytes())
#free all remaining variables (phi_k_pad, k_arr, kernel_output, xx_flat, zz_flat, yy_flat, out_records, kernel_real, kernel_imag)

