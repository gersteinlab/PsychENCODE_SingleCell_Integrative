��
��
D
AddV2
x"T
y"T
z"T"
Ttype:
2	��
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( �
�
BiasAdd

value"T	
bias"T
output"T""
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
$
DisableCopyOnRead
resource�
.
Identity

input"T
output"T"	
Ttype
u
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:
2	
�
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool("
allow_missing_filesbool( �
?
Mul
x"T
y"T
z"T"
Ttype:
2	�

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
.
Rsqrt
x"T
y"T"
Ttype:

2
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
0
Sigmoid
x"T
y"T"
Ttype:

2
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ��
@
StaticRegexFullMatch	
input

output
"
patternstring
L

StringJoin
inputs*N

output"

Nint("
	separatorstring 
<
Sub
x"T
y"T
z"T"
Ttype:
2	
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.13.02v2.13.0-rc2-7-g1cb1a030a628��

^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
�
$Adam/v/dense_model_42/dense_128/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/v/dense_model_42/dense_128/bias
�
8Adam/v/dense_model_42/dense_128/bias/Read/ReadVariableOpReadVariableOp$Adam/v/dense_model_42/dense_128/bias*
_output_shapes
:*
dtype0
�
$Adam/m/dense_model_42/dense_128/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/m/dense_model_42/dense_128/bias
�
8Adam/m/dense_model_42/dense_128/bias/Read/ReadVariableOpReadVariableOp$Adam/m/dense_model_42/dense_128/bias*
_output_shapes
:*
dtype0
�
&Adam/v/dense_model_42/dense_128/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/v/dense_model_42/dense_128/kernel
�
:Adam/v/dense_model_42/dense_128/kernel/Read/ReadVariableOpReadVariableOp&Adam/v/dense_model_42/dense_128/kernel*
_output_shapes

:*
dtype0
�
&Adam/m/dense_model_42/dense_128/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/m/dense_model_42/dense_128/kernel
�
:Adam/m/dense_model_42/dense_128/kernel/Read/ReadVariableOpReadVariableOp&Adam/m/dense_model_42/dense_128/kernel*
_output_shapes

:*
dtype0
�
1Adam/v/dense_model_42/batch_normalization_85/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31Adam/v/dense_model_42/batch_normalization_85/beta
�
EAdam/v/dense_model_42/batch_normalization_85/beta/Read/ReadVariableOpReadVariableOp1Adam/v/dense_model_42/batch_normalization_85/beta*
_output_shapes
:*
dtype0
�
1Adam/m/dense_model_42/batch_normalization_85/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31Adam/m/dense_model_42/batch_normalization_85/beta
�
EAdam/m/dense_model_42/batch_normalization_85/beta/Read/ReadVariableOpReadVariableOp1Adam/m/dense_model_42/batch_normalization_85/beta*
_output_shapes
:*
dtype0
�
2Adam/v/dense_model_42/batch_normalization_85/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/v/dense_model_42/batch_normalization_85/gamma
�
FAdam/v/dense_model_42/batch_normalization_85/gamma/Read/ReadVariableOpReadVariableOp2Adam/v/dense_model_42/batch_normalization_85/gamma*
_output_shapes
:*
dtype0
�
2Adam/m/dense_model_42/batch_normalization_85/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/m/dense_model_42/batch_normalization_85/gamma
�
FAdam/m/dense_model_42/batch_normalization_85/gamma/Read/ReadVariableOpReadVariableOp2Adam/m/dense_model_42/batch_normalization_85/gamma*
_output_shapes
:*
dtype0
�
$Adam/v/dense_model_42/dense_127/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/v/dense_model_42/dense_127/bias
�
8Adam/v/dense_model_42/dense_127/bias/Read/ReadVariableOpReadVariableOp$Adam/v/dense_model_42/dense_127/bias*
_output_shapes
:*
dtype0
�
$Adam/m/dense_model_42/dense_127/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/m/dense_model_42/dense_127/bias
�
8Adam/m/dense_model_42/dense_127/bias/Read/ReadVariableOpReadVariableOp$Adam/m/dense_model_42/dense_127/bias*
_output_shapes
:*
dtype0
�
&Adam/v/dense_model_42/dense_127/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/v/dense_model_42/dense_127/kernel
�
:Adam/v/dense_model_42/dense_127/kernel/Read/ReadVariableOpReadVariableOp&Adam/v/dense_model_42/dense_127/kernel*
_output_shapes

:*
dtype0
�
&Adam/m/dense_model_42/dense_127/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/m/dense_model_42/dense_127/kernel
�
:Adam/m/dense_model_42/dense_127/kernel/Read/ReadVariableOpReadVariableOp&Adam/m/dense_model_42/dense_127/kernel*
_output_shapes

:*
dtype0
�
1Adam/v/dense_model_42/batch_normalization_84/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31Adam/v/dense_model_42/batch_normalization_84/beta
�
EAdam/v/dense_model_42/batch_normalization_84/beta/Read/ReadVariableOpReadVariableOp1Adam/v/dense_model_42/batch_normalization_84/beta*
_output_shapes
:*
dtype0
�
1Adam/m/dense_model_42/batch_normalization_84/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31Adam/m/dense_model_42/batch_normalization_84/beta
�
EAdam/m/dense_model_42/batch_normalization_84/beta/Read/ReadVariableOpReadVariableOp1Adam/m/dense_model_42/batch_normalization_84/beta*
_output_shapes
:*
dtype0
�
2Adam/v/dense_model_42/batch_normalization_84/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/v/dense_model_42/batch_normalization_84/gamma
�
FAdam/v/dense_model_42/batch_normalization_84/gamma/Read/ReadVariableOpReadVariableOp2Adam/v/dense_model_42/batch_normalization_84/gamma*
_output_shapes
:*
dtype0
�
2Adam/m/dense_model_42/batch_normalization_84/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/m/dense_model_42/batch_normalization_84/gamma
�
FAdam/m/dense_model_42/batch_normalization_84/gamma/Read/ReadVariableOpReadVariableOp2Adam/m/dense_model_42/batch_normalization_84/gamma*
_output_shapes
:*
dtype0
�
$Adam/v/dense_model_42/dense_126/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/v/dense_model_42/dense_126/bias
�
8Adam/v/dense_model_42/dense_126/bias/Read/ReadVariableOpReadVariableOp$Adam/v/dense_model_42/dense_126/bias*
_output_shapes
:*
dtype0
�
$Adam/m/dense_model_42/dense_126/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/m/dense_model_42/dense_126/bias
�
8Adam/m/dense_model_42/dense_126/bias/Read/ReadVariableOpReadVariableOp$Adam/m/dense_model_42/dense_126/bias*
_output_shapes
:*
dtype0
�
&Adam/v/dense_model_42/dense_126/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*7
shared_name(&Adam/v/dense_model_42/dense_126/kernel
�
:Adam/v/dense_model_42/dense_126/kernel/Read/ReadVariableOpReadVariableOp&Adam/v/dense_model_42/dense_126/kernel*
_output_shapes
:	�*
dtype0
�
&Adam/m/dense_model_42/dense_126/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*7
shared_name(&Adam/m/dense_model_42/dense_126/kernel
�
:Adam/m/dense_model_42/dense_126/kernel/Read/ReadVariableOpReadVariableOp&Adam/m/dense_model_42/dense_126/kernel*
_output_shapes
:	�*
dtype0
n
learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namelearning_rate
g
!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate*
_output_shapes
: *
dtype0
f
	iterationVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	iteration
_
iteration/Read/ReadVariableOpReadVariableOp	iteration*
_output_shapes
: *
dtype0	
�
dense_model_42/dense_128/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*.
shared_namedense_model_42/dense_128/bias
�
1dense_model_42/dense_128/bias/Read/ReadVariableOpReadVariableOpdense_model_42/dense_128/bias*
_output_shapes
:*
dtype0
�
dense_model_42/dense_128/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*0
shared_name!dense_model_42/dense_128/kernel
�
3dense_model_42/dense_128/kernel/Read/ReadVariableOpReadVariableOpdense_model_42/dense_128/kernel*
_output_shapes

:*
dtype0
�
5dense_model_42/batch_normalization_85/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*F
shared_name75dense_model_42/batch_normalization_85/moving_variance
�
Idense_model_42/batch_normalization_85/moving_variance/Read/ReadVariableOpReadVariableOp5dense_model_42/batch_normalization_85/moving_variance*
_output_shapes
:*
dtype0
�
1dense_model_42/batch_normalization_85/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31dense_model_42/batch_normalization_85/moving_mean
�
Edense_model_42/batch_normalization_85/moving_mean/Read/ReadVariableOpReadVariableOp1dense_model_42/batch_normalization_85/moving_mean*
_output_shapes
:*
dtype0
�
*dense_model_42/batch_normalization_85/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*;
shared_name,*dense_model_42/batch_normalization_85/beta
�
>dense_model_42/batch_normalization_85/beta/Read/ReadVariableOpReadVariableOp*dense_model_42/batch_normalization_85/beta*
_output_shapes
:*
dtype0
�
+dense_model_42/batch_normalization_85/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*<
shared_name-+dense_model_42/batch_normalization_85/gamma
�
?dense_model_42/batch_normalization_85/gamma/Read/ReadVariableOpReadVariableOp+dense_model_42/batch_normalization_85/gamma*
_output_shapes
:*
dtype0
�
dense_model_42/dense_127/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*.
shared_namedense_model_42/dense_127/bias
�
1dense_model_42/dense_127/bias/Read/ReadVariableOpReadVariableOpdense_model_42/dense_127/bias*
_output_shapes
:*
dtype0
�
dense_model_42/dense_127/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*0
shared_name!dense_model_42/dense_127/kernel
�
3dense_model_42/dense_127/kernel/Read/ReadVariableOpReadVariableOpdense_model_42/dense_127/kernel*
_output_shapes

:*
dtype0
�
5dense_model_42/batch_normalization_84/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*F
shared_name75dense_model_42/batch_normalization_84/moving_variance
�
Idense_model_42/batch_normalization_84/moving_variance/Read/ReadVariableOpReadVariableOp5dense_model_42/batch_normalization_84/moving_variance*
_output_shapes
:*
dtype0
�
1dense_model_42/batch_normalization_84/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31dense_model_42/batch_normalization_84/moving_mean
�
Edense_model_42/batch_normalization_84/moving_mean/Read/ReadVariableOpReadVariableOp1dense_model_42/batch_normalization_84/moving_mean*
_output_shapes
:*
dtype0
�
*dense_model_42/batch_normalization_84/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*;
shared_name,*dense_model_42/batch_normalization_84/beta
�
>dense_model_42/batch_normalization_84/beta/Read/ReadVariableOpReadVariableOp*dense_model_42/batch_normalization_84/beta*
_output_shapes
:*
dtype0
�
+dense_model_42/batch_normalization_84/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*<
shared_name-+dense_model_42/batch_normalization_84/gamma
�
?dense_model_42/batch_normalization_84/gamma/Read/ReadVariableOpReadVariableOp+dense_model_42/batch_normalization_84/gamma*
_output_shapes
:*
dtype0
�
dense_model_42/dense_126/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*.
shared_namedense_model_42/dense_126/bias
�
1dense_model_42/dense_126/bias/Read/ReadVariableOpReadVariableOpdense_model_42/dense_126/bias*
_output_shapes
:*
dtype0
�
dense_model_42/dense_126/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*0
shared_name!dense_model_42/dense_126/kernel
�
3dense_model_42/dense_126/kernel/Read/ReadVariableOpReadVariableOpdense_model_42/dense_126/kernel*
_output_shapes
:	�*
dtype0
|
serving_default_input_1Placeholder*(
_output_shapes
:����������*
dtype0*
shape:����������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1dense_model_42/dense_126/kerneldense_model_42/dense_126/bias1dense_model_42/batch_normalization_84/moving_mean5dense_model_42/batch_normalization_84/moving_variance*dense_model_42/batch_normalization_84/beta+dense_model_42/batch_normalization_84/gammadense_model_42/dense_127/kerneldense_model_42/dense_127/bias1dense_model_42/batch_normalization_85/moving_mean5dense_model_42/batch_normalization_85/moving_variance*dense_model_42/batch_normalization_85/beta+dense_model_42/batch_normalization_85/gammadense_model_42/dense_128/kerneldense_model_42/dense_128/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *-
f(R&
$__inference_signature_wrapper_161720

NoOpNoOp
�G
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�G
value�FB�F B�F
�
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature

dense1
		norm1


dense2
	norm2
out
	optimizer

signatures*
j
0
1
2
3
4
5
6
7
8
9
10
11
12
13*
J
0
1
2
3
4
5
6
7
8
9*
* 
�
non_trainable_variables

layers
metrics
 layer_regularization_losses
!layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*

"trace_0
#trace_1* 

$trace_0
%trace_1* 
* 
�
&	variables
'trainable_variables
(regularization_losses
)	keras_api
*__call__
*+&call_and_return_all_conditional_losses

kernel
bias*
�
,	variables
-trainable_variables
.regularization_losses
/	keras_api
0__call__
*1&call_and_return_all_conditional_losses
2axis
	gamma
beta
moving_mean
moving_variance*
�
3	variables
4trainable_variables
5regularization_losses
6	keras_api
7__call__
*8&call_and_return_all_conditional_losses

kernel
bias*
�
9	variables
:trainable_variables
;regularization_losses
<	keras_api
=__call__
*>&call_and_return_all_conditional_losses
?axis
	gamma
beta
moving_mean
moving_variance*
�
@	variables
Atrainable_variables
Bregularization_losses
C	keras_api
D__call__
*E&call_and_return_all_conditional_losses

kernel
bias*
�
F
_variables
G_iterations
H_learning_rate
I_index_dict
J
_momentums
K_velocities
L_update_step_xla*

Mserving_default* 
_Y
VARIABLE_VALUEdense_model_42/dense_126/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEdense_model_42/dense_126/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
ke
VARIABLE_VALUE+dense_model_42/batch_normalization_84/gamma&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
jd
VARIABLE_VALUE*dense_model_42/batch_normalization_84/beta&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE1dense_model_42/batch_normalization_84/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
uo
VARIABLE_VALUE5dense_model_42/batch_normalization_84/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEdense_model_42/dense_127/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEdense_model_42/dense_127/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE*
ke
VARIABLE_VALUE+dense_model_42/batch_normalization_85/gamma&variables/8/.ATTRIBUTES/VARIABLE_VALUE*
jd
VARIABLE_VALUE*dense_model_42/batch_normalization_85/beta&variables/9/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE1dense_model_42/batch_normalization_85/moving_mean'variables/10/.ATTRIBUTES/VARIABLE_VALUE*
vp
VARIABLE_VALUE5dense_model_42/batch_normalization_85/moving_variance'variables/11/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEdense_model_42/dense_128/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEdense_model_42/dense_128/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE*
 
0
1
2
3*
'
0
	1

2
3
4*

N0
O1*
* 
* 
* 
* 
* 
* 

0
1*

0
1*
* 
�
Pnon_trainable_variables

Qlayers
Rmetrics
Slayer_regularization_losses
Tlayer_metrics
&	variables
'trainable_variables
(regularization_losses
*__call__
*+&call_and_return_all_conditional_losses
&+"call_and_return_conditional_losses*

Utrace_0* 

Vtrace_0* 
 
0
1
2
3*

0
1*
* 
�
Wnon_trainable_variables

Xlayers
Ymetrics
Zlayer_regularization_losses
[layer_metrics
,	variables
-trainable_variables
.regularization_losses
0__call__
*1&call_and_return_all_conditional_losses
&1"call_and_return_conditional_losses*

\trace_0
]trace_1* 

^trace_0
_trace_1* 
* 

0
1*

0
1*
* 
�
`non_trainable_variables

alayers
bmetrics
clayer_regularization_losses
dlayer_metrics
3	variables
4trainable_variables
5regularization_losses
7__call__
*8&call_and_return_all_conditional_losses
&8"call_and_return_conditional_losses*

etrace_0* 

ftrace_0* 
 
0
1
2
3*

0
1*
* 
�
gnon_trainable_variables

hlayers
imetrics
jlayer_regularization_losses
klayer_metrics
9	variables
:trainable_variables
;regularization_losses
=__call__
*>&call_and_return_all_conditional_losses
&>"call_and_return_conditional_losses*

ltrace_0
mtrace_1* 

ntrace_0
otrace_1* 
* 

0
1*

0
1*
* 
�
pnon_trainable_variables

qlayers
rmetrics
slayer_regularization_losses
tlayer_metrics
@	variables
Atrainable_variables
Bregularization_losses
D__call__
*E&call_and_return_all_conditional_losses
&E"call_and_return_conditional_losses*

utrace_0* 

vtrace_0* 
�
G0
w1
x2
y3
z4
{5
|6
}7
~8
9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20*
SM
VARIABLE_VALUE	iteration0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUElearning_rate3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 
O
w0
y1
{2
}3
4
�5
�6
�7
�8
�9*
P
x0
z1
|2
~3
�4
�5
�6
�7
�8
�9*
* 
* 
<
�	variables
�	keras_api

�total

�count*
M
�	variables
�	keras_api

�total

�count
�
_fn_kwargs*
* 
* 
* 
* 
* 
* 
* 

0
1*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 

0
1*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
qk
VARIABLE_VALUE&Adam/m/dense_model_42/dense_126/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE&Adam/v/dense_model_42/dense_126/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
oi
VARIABLE_VALUE$Adam/m/dense_model_42/dense_126/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
oi
VARIABLE_VALUE$Adam/v/dense_model_42/dense_126/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
}w
VARIABLE_VALUE2Adam/m/dense_model_42/batch_normalization_84/gamma1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
}w
VARIABLE_VALUE2Adam/v/dense_model_42/batch_normalization_84/gamma1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
|v
VARIABLE_VALUE1Adam/m/dense_model_42/batch_normalization_84/beta1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
|v
VARIABLE_VALUE1Adam/v/dense_model_42/batch_normalization_84/beta1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE&Adam/m/dense_model_42/dense_127/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/v/dense_model_42/dense_127/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/m/dense_model_42/dense_127/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/v/dense_model_42/dense_127/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE2Adam/m/dense_model_42/batch_normalization_85/gamma2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE2Adam/v/dense_model_42/batch_normalization_85/gamma2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUE*
}w
VARIABLE_VALUE1Adam/m/dense_model_42/batch_normalization_85/beta2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUE*
}w
VARIABLE_VALUE1Adam/v/dense_model_42/batch_normalization_85/beta2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/m/dense_model_42/dense_128/kernel2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/v/dense_model_42/dense_128/kernel2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/m/dense_model_42/dense_128/bias2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/v/dense_model_42/dense_128/bias2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
SM
VARIABLE_VALUEtotal4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE*
* 
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenamedense_model_42/dense_126/kerneldense_model_42/dense_126/bias+dense_model_42/batch_normalization_84/gamma*dense_model_42/batch_normalization_84/beta1dense_model_42/batch_normalization_84/moving_mean5dense_model_42/batch_normalization_84/moving_variancedense_model_42/dense_127/kerneldense_model_42/dense_127/bias+dense_model_42/batch_normalization_85/gamma*dense_model_42/batch_normalization_85/beta1dense_model_42/batch_normalization_85/moving_mean5dense_model_42/batch_normalization_85/moving_variancedense_model_42/dense_128/kerneldense_model_42/dense_128/bias	iterationlearning_rate&Adam/m/dense_model_42/dense_126/kernel&Adam/v/dense_model_42/dense_126/kernel$Adam/m/dense_model_42/dense_126/bias$Adam/v/dense_model_42/dense_126/bias2Adam/m/dense_model_42/batch_normalization_84/gamma2Adam/v/dense_model_42/batch_normalization_84/gamma1Adam/m/dense_model_42/batch_normalization_84/beta1Adam/v/dense_model_42/batch_normalization_84/beta&Adam/m/dense_model_42/dense_127/kernel&Adam/v/dense_model_42/dense_127/kernel$Adam/m/dense_model_42/dense_127/bias$Adam/v/dense_model_42/dense_127/bias2Adam/m/dense_model_42/batch_normalization_85/gamma2Adam/v/dense_model_42/batch_normalization_85/gamma1Adam/m/dense_model_42/batch_normalization_85/beta1Adam/v/dense_model_42/batch_normalization_85/beta&Adam/m/dense_model_42/dense_128/kernel&Adam/v/dense_model_42/dense_128/kernel$Adam/m/dense_model_42/dense_128/bias$Adam/v/dense_model_42/dense_128/biastotal_1count_1totalcountConst*5
Tin.
,2**
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *(
f#R!
__inference__traced_save_162200
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_model_42/dense_126/kerneldense_model_42/dense_126/bias+dense_model_42/batch_normalization_84/gamma*dense_model_42/batch_normalization_84/beta1dense_model_42/batch_normalization_84/moving_mean5dense_model_42/batch_normalization_84/moving_variancedense_model_42/dense_127/kerneldense_model_42/dense_127/bias+dense_model_42/batch_normalization_85/gamma*dense_model_42/batch_normalization_85/beta1dense_model_42/batch_normalization_85/moving_mean5dense_model_42/batch_normalization_85/moving_variancedense_model_42/dense_128/kerneldense_model_42/dense_128/bias	iterationlearning_rate&Adam/m/dense_model_42/dense_126/kernel&Adam/v/dense_model_42/dense_126/kernel$Adam/m/dense_model_42/dense_126/bias$Adam/v/dense_model_42/dense_126/bias2Adam/m/dense_model_42/batch_normalization_84/gamma2Adam/v/dense_model_42/batch_normalization_84/gamma1Adam/m/dense_model_42/batch_normalization_84/beta1Adam/v/dense_model_42/batch_normalization_84/beta&Adam/m/dense_model_42/dense_127/kernel&Adam/v/dense_model_42/dense_127/kernel$Adam/m/dense_model_42/dense_127/bias$Adam/v/dense_model_42/dense_127/bias2Adam/m/dense_model_42/batch_normalization_85/gamma2Adam/v/dense_model_42/batch_normalization_85/gamma1Adam/m/dense_model_42/batch_normalization_85/beta1Adam/v/dense_model_42/batch_normalization_85/beta&Adam/m/dense_model_42/dense_128/kernel&Adam/v/dense_model_42/dense_128/kernel$Adam/m/dense_model_42/dense_128/bias$Adam/v/dense_model_42/dense_128/biastotal_1count_1totalcount*4
Tin-
+2)*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *+
f&R$
"__inference__traced_restore_162329��
��
�
"__inference__traced_restore_162329
file_prefixC
0assignvariableop_dense_model_42_dense_126_kernel:	�>
0assignvariableop_1_dense_model_42_dense_126_bias:L
>assignvariableop_2_dense_model_42_batch_normalization_84_gamma:K
=assignvariableop_3_dense_model_42_batch_normalization_84_beta:R
Dassignvariableop_4_dense_model_42_batch_normalization_84_moving_mean:V
Hassignvariableop_5_dense_model_42_batch_normalization_84_moving_variance:D
2assignvariableop_6_dense_model_42_dense_127_kernel:>
0assignvariableop_7_dense_model_42_dense_127_bias:L
>assignvariableop_8_dense_model_42_batch_normalization_85_gamma:K
=assignvariableop_9_dense_model_42_batch_normalization_85_beta:S
Eassignvariableop_10_dense_model_42_batch_normalization_85_moving_mean:W
Iassignvariableop_11_dense_model_42_batch_normalization_85_moving_variance:E
3assignvariableop_12_dense_model_42_dense_128_kernel:?
1assignvariableop_13_dense_model_42_dense_128_bias:'
assignvariableop_14_iteration:	 +
!assignvariableop_15_learning_rate: M
:assignvariableop_16_adam_m_dense_model_42_dense_126_kernel:	�M
:assignvariableop_17_adam_v_dense_model_42_dense_126_kernel:	�F
8assignvariableop_18_adam_m_dense_model_42_dense_126_bias:F
8assignvariableop_19_adam_v_dense_model_42_dense_126_bias:T
Fassignvariableop_20_adam_m_dense_model_42_batch_normalization_84_gamma:T
Fassignvariableop_21_adam_v_dense_model_42_batch_normalization_84_gamma:S
Eassignvariableop_22_adam_m_dense_model_42_batch_normalization_84_beta:S
Eassignvariableop_23_adam_v_dense_model_42_batch_normalization_84_beta:L
:assignvariableop_24_adam_m_dense_model_42_dense_127_kernel:L
:assignvariableop_25_adam_v_dense_model_42_dense_127_kernel:F
8assignvariableop_26_adam_m_dense_model_42_dense_127_bias:F
8assignvariableop_27_adam_v_dense_model_42_dense_127_bias:T
Fassignvariableop_28_adam_m_dense_model_42_batch_normalization_85_gamma:T
Fassignvariableop_29_adam_v_dense_model_42_batch_normalization_85_gamma:S
Eassignvariableop_30_adam_m_dense_model_42_batch_normalization_85_beta:S
Eassignvariableop_31_adam_v_dense_model_42_batch_normalization_85_beta:L
:assignvariableop_32_adam_m_dense_model_42_dense_128_kernel:L
:assignvariableop_33_adam_v_dense_model_42_dense_128_kernel:F
8assignvariableop_34_adam_m_dense_model_42_dense_128_bias:F
8assignvariableop_35_adam_v_dense_model_42_dense_128_bias:%
assignvariableop_36_total_1: %
assignvariableop_37_count_1: #
assignvariableop_38_total: #
assignvariableop_39_count: 
identity_41��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:)*
dtype0*�
value�B�)B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:)*
dtype0*e
value\BZ)B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�:::::::::::::::::::::::::::::::::::::::::*7
dtypes-
+2)	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOp0assignvariableop_dense_model_42_dense_126_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp0assignvariableop_1_dense_model_42_dense_126_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOp>assignvariableop_2_dense_model_42_batch_normalization_84_gammaIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOp=assignvariableop_3_dense_model_42_batch_normalization_84_betaIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpDassignvariableop_4_dense_model_42_batch_normalization_84_moving_meanIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpHassignvariableop_5_dense_model_42_batch_normalization_84_moving_varianceIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp2assignvariableop_6_dense_model_42_dense_127_kernelIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp0assignvariableop_7_dense_model_42_dense_127_biasIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOp>assignvariableop_8_dense_model_42_batch_normalization_85_gammaIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp=assignvariableop_9_dense_model_42_batch_normalization_85_betaIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOpEassignvariableop_10_dense_model_42_batch_normalization_85_moving_meanIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpIassignvariableop_11_dense_model_42_batch_normalization_85_moving_varianceIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp3assignvariableop_12_dense_model_42_dense_128_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp1assignvariableop_13_dense_model_42_dense_128_biasIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_14AssignVariableOpassignvariableop_14_iterationIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOp!assignvariableop_15_learning_rateIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOp:assignvariableop_16_adam_m_dense_model_42_dense_126_kernelIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOp:assignvariableop_17_adam_v_dense_model_42_dense_126_kernelIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp8assignvariableop_18_adam_m_dense_model_42_dense_126_biasIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp8assignvariableop_19_adam_v_dense_model_42_dense_126_biasIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOpFassignvariableop_20_adam_m_dense_model_42_batch_normalization_84_gammaIdentity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOpFassignvariableop_21_adam_v_dense_model_42_batch_normalization_84_gammaIdentity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOpEassignvariableop_22_adam_m_dense_model_42_batch_normalization_84_betaIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOpEassignvariableop_23_adam_v_dense_model_42_batch_normalization_84_betaIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOp:assignvariableop_24_adam_m_dense_model_42_dense_127_kernelIdentity_24:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOp:assignvariableop_25_adam_v_dense_model_42_dense_127_kernelIdentity_25:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOp8assignvariableop_26_adam_m_dense_model_42_dense_127_biasIdentity_26:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOp8assignvariableop_27_adam_v_dense_model_42_dense_127_biasIdentity_27:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpFassignvariableop_28_adam_m_dense_model_42_batch_normalization_85_gammaIdentity_28:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpFassignvariableop_29_adam_v_dense_model_42_batch_normalization_85_gammaIdentity_29:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOpEassignvariableop_30_adam_m_dense_model_42_batch_normalization_85_betaIdentity_30:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOpEassignvariableop_31_adam_v_dense_model_42_batch_normalization_85_betaIdentity_31:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOp:assignvariableop_32_adam_m_dense_model_42_dense_128_kernelIdentity_32:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOp:assignvariableop_33_adam_v_dense_model_42_dense_128_kernelIdentity_33:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOp8assignvariableop_34_adam_m_dense_model_42_dense_128_biasIdentity_34:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOp8assignvariableop_35_adam_v_dense_model_42_dense_128_biasIdentity_35:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_36AssignVariableOpassignvariableop_36_total_1Identity_36:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_37AssignVariableOpassignvariableop_37_count_1Identity_37:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_38AssignVariableOpassignvariableop_38_totalIdentity_38:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_39AssignVariableOpassignvariableop_39_countIdentity_39:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0Y
NoOpNoOp"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 �
Identity_40Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_41IdentityIdentity_40:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
_output_shapes
 "#
identity_41Identity_41:output:0*(
_construction_contextkEagerRuntime*e
_input_shapesT
R: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92$
AssignVariableOpAssignVariableOp:%(!

_user_specified_namecount:%'!

_user_specified_nametotal:'&#
!
_user_specified_name	count_1:'%#
!
_user_specified_name	total_1:D$@
>
_user_specified_name&$Adam/v/dense_model_42/dense_128/bias:D#@
>
_user_specified_name&$Adam/m/dense_model_42/dense_128/bias:F"B
@
_user_specified_name(&Adam/v/dense_model_42/dense_128/kernel:F!B
@
_user_specified_name(&Adam/m/dense_model_42/dense_128/kernel:Q M
K
_user_specified_name31Adam/v/dense_model_42/batch_normalization_85/beta:QM
K
_user_specified_name31Adam/m/dense_model_42/batch_normalization_85/beta:RN
L
_user_specified_name42Adam/v/dense_model_42/batch_normalization_85/gamma:RN
L
_user_specified_name42Adam/m/dense_model_42/batch_normalization_85/gamma:D@
>
_user_specified_name&$Adam/v/dense_model_42/dense_127/bias:D@
>
_user_specified_name&$Adam/m/dense_model_42/dense_127/bias:FB
@
_user_specified_name(&Adam/v/dense_model_42/dense_127/kernel:FB
@
_user_specified_name(&Adam/m/dense_model_42/dense_127/kernel:QM
K
_user_specified_name31Adam/v/dense_model_42/batch_normalization_84/beta:QM
K
_user_specified_name31Adam/m/dense_model_42/batch_normalization_84/beta:RN
L
_user_specified_name42Adam/v/dense_model_42/batch_normalization_84/gamma:RN
L
_user_specified_name42Adam/m/dense_model_42/batch_normalization_84/gamma:D@
>
_user_specified_name&$Adam/v/dense_model_42/dense_126/bias:D@
>
_user_specified_name&$Adam/m/dense_model_42/dense_126/bias:FB
@
_user_specified_name(&Adam/v/dense_model_42/dense_126/kernel:FB
@
_user_specified_name(&Adam/m/dense_model_42/dense_126/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:=9
7
_user_specified_namedense_model_42/dense_128/bias:?;
9
_user_specified_name!dense_model_42/dense_128/kernel:UQ
O
_user_specified_name75dense_model_42/batch_normalization_85/moving_variance:QM
K
_user_specified_name31dense_model_42/batch_normalization_85/moving_mean:J
F
D
_user_specified_name,*dense_model_42/batch_normalization_85/beta:K	G
E
_user_specified_name-+dense_model_42/batch_normalization_85/gamma:=9
7
_user_specified_namedense_model_42/dense_127/bias:?;
9
_user_specified_name!dense_model_42/dense_127/kernel:UQ
O
_user_specified_name75dense_model_42/batch_normalization_84/moving_variance:QM
K
_user_specified_name31dense_model_42/batch_normalization_84/moving_mean:JF
D
_user_specified_name,*dense_model_42/batch_normalization_84/beta:KG
E
_user_specified_name-+dense_model_42/batch_normalization_84/gamma:=9
7
_user_specified_namedense_model_42/dense_126/bias:?;
9
_user_specified_name!dense_model_42/dense_126/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161460

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:t
batchnorm/addAddV2Cast_1/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�$
�
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161591
input_1#
dense_126_161557:	�
dense_126_161559:+
batch_normalization_84_161562:+
batch_normalization_84_161564:+
batch_normalization_84_161566:+
batch_normalization_84_161568:"
dense_127_161571:
dense_127_161573:+
batch_normalization_85_161576:+
batch_normalization_85_161578:+
batch_normalization_85_161580:+
batch_normalization_85_161582:"
dense_128_161585:
dense_128_161587:
identity��.batch_normalization_84/StatefulPartitionedCall�.batch_normalization_85/StatefulPartitionedCall�!dense_126/StatefulPartitionedCall�!dense_127/StatefulPartitionedCall�!dense_128/StatefulPartitionedCall�
!dense_126/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_126_161557dense_126_161559*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_126_layer_call_and_return_conditional_losses_161498�
.batch_normalization_84/StatefulPartitionedCallStatefulPartitionedCall*dense_126/StatefulPartitionedCall:output:0batch_normalization_84_161562batch_normalization_84_161564batch_normalization_84_161566batch_normalization_84_161568*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161380�
!dense_127/StatefulPartitionedCallStatefulPartitionedCall7batch_normalization_84/StatefulPartitionedCall:output:0dense_127_161571dense_127_161573*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_127_layer_call_and_return_conditional_losses_161522�
.batch_normalization_85/StatefulPartitionedCallStatefulPartitionedCall*dense_127/StatefulPartitionedCall:output:0batch_normalization_85_161576batch_normalization_85_161578batch_normalization_85_161580batch_normalization_85_161582*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161460�
!dense_128/StatefulPartitionedCallStatefulPartitionedCall7batch_normalization_85/StatefulPartitionedCall:output:0dense_128_161585dense_128_161587*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_128_layer_call_and_return_conditional_losses_161547y
IdentityIdentity*dense_128/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp/^batch_normalization_84/StatefulPartitionedCall/^batch_normalization_85/StatefulPartitionedCall"^dense_126/StatefulPartitionedCall"^dense_127/StatefulPartitionedCall"^dense_128/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2`
.batch_normalization_84/StatefulPartitionedCall.batch_normalization_84/StatefulPartitionedCall2`
.batch_normalization_85/StatefulPartitionedCall.batch_normalization_85/StatefulPartitionedCall2F
!dense_126/StatefulPartitionedCall!dense_126/StatefulPartitionedCall2F
!dense_127/StatefulPartitionedCall!dense_127/StatefulPartitionedCall2F
!dense_128/StatefulPartitionedCall!dense_128/StatefulPartitionedCall:&"
 
_user_specified_name161587:&"
 
_user_specified_name161585:&"
 
_user_specified_name161582:&"
 
_user_specified_name161580:&
"
 
_user_specified_name161578:&	"
 
_user_specified_name161576:&"
 
_user_specified_name161573:&"
 
_user_specified_name161571:&"
 
_user_specified_name161568:&"
 
_user_specified_name161566:&"
 
_user_specified_name161564:&"
 
_user_specified_name161562:&"
 
_user_specified_name161559:&"
 
_user_specified_name161557:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161918

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:t
batchnorm/addAddV2Cast_1/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
7__inference_batch_normalization_85_layer_call_fn_161851

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161440o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161847:&"
 
_user_specified_name161845:&"
 
_user_specified_name161843:&"
 
_user_specified_name161841:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_126_layer_call_and_return_conditional_losses_161498

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
/__inference_dense_model_42_layer_call_fn_161624
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*,
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161554o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161620:&"
 
_user_specified_name161618:&"
 
_user_specified_name161616:&"
 
_user_specified_name161614:&
"
 
_user_specified_name161612:&	"
 
_user_specified_name161610:&"
 
_user_specified_name161608:&"
 
_user_specified_name161606:&"
 
_user_specified_name161604:&"
 
_user_specified_name161602:&"
 
_user_specified_name161600:&"
 
_user_specified_name161598:&"
 
_user_specified_name161596:&"
 
_user_specified_name161594:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
*__inference_dense_126_layer_call_fn_161729

inputs
unknown:	�
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_126_layer_call_and_return_conditional_losses_161498o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161725:&"
 
_user_specified_name161723:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�

�
E__inference_dense_128_layer_call_and_return_conditional_losses_161938

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:���������Z
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161380

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:t
batchnorm/addAddV2Cast_1/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�%
�
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161799

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 Z
AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
AssignMovingAvg/ReadVariableOpReadVariableOp'assignmovingavg_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvgAssignSubVariableOp'assignmovingavg_readvariableop_resourceAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0\
AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
 AssignMovingAvg_1/ReadVariableOpReadVariableOp)assignmovingavg_1_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:q
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12"
AssignMovingAvgAssignMovingAvg2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�%
�
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161360

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 Z
AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
AssignMovingAvg/ReadVariableOpReadVariableOp'assignmovingavg_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvgAssignSubVariableOp'assignmovingavg_readvariableop_resourceAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0\
AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
 AssignMovingAvg_1/ReadVariableOpReadVariableOp)assignmovingavg_1_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:q
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12"
AssignMovingAvgAssignMovingAvg2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_127_layer_call_and_return_conditional_losses_161522

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
7__inference_batch_normalization_84_layer_call_fn_161765

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161380o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161761:&"
 
_user_specified_name161759:&"
 
_user_specified_name161757:&"
 
_user_specified_name161755:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
��
�)
__inference__traced_save_162200
file_prefixI
6read_disablecopyonread_dense_model_42_dense_126_kernel:	�D
6read_1_disablecopyonread_dense_model_42_dense_126_bias:R
Dread_2_disablecopyonread_dense_model_42_batch_normalization_84_gamma:Q
Cread_3_disablecopyonread_dense_model_42_batch_normalization_84_beta:X
Jread_4_disablecopyonread_dense_model_42_batch_normalization_84_moving_mean:\
Nread_5_disablecopyonread_dense_model_42_batch_normalization_84_moving_variance:J
8read_6_disablecopyonread_dense_model_42_dense_127_kernel:D
6read_7_disablecopyonread_dense_model_42_dense_127_bias:R
Dread_8_disablecopyonread_dense_model_42_batch_normalization_85_gamma:Q
Cread_9_disablecopyonread_dense_model_42_batch_normalization_85_beta:Y
Kread_10_disablecopyonread_dense_model_42_batch_normalization_85_moving_mean:]
Oread_11_disablecopyonread_dense_model_42_batch_normalization_85_moving_variance:K
9read_12_disablecopyonread_dense_model_42_dense_128_kernel:E
7read_13_disablecopyonread_dense_model_42_dense_128_bias:-
#read_14_disablecopyonread_iteration:	 1
'read_15_disablecopyonread_learning_rate: S
@read_16_disablecopyonread_adam_m_dense_model_42_dense_126_kernel:	�S
@read_17_disablecopyonread_adam_v_dense_model_42_dense_126_kernel:	�L
>read_18_disablecopyonread_adam_m_dense_model_42_dense_126_bias:L
>read_19_disablecopyonread_adam_v_dense_model_42_dense_126_bias:Z
Lread_20_disablecopyonread_adam_m_dense_model_42_batch_normalization_84_gamma:Z
Lread_21_disablecopyonread_adam_v_dense_model_42_batch_normalization_84_gamma:Y
Kread_22_disablecopyonread_adam_m_dense_model_42_batch_normalization_84_beta:Y
Kread_23_disablecopyonread_adam_v_dense_model_42_batch_normalization_84_beta:R
@read_24_disablecopyonread_adam_m_dense_model_42_dense_127_kernel:R
@read_25_disablecopyonread_adam_v_dense_model_42_dense_127_kernel:L
>read_26_disablecopyonread_adam_m_dense_model_42_dense_127_bias:L
>read_27_disablecopyonread_adam_v_dense_model_42_dense_127_bias:Z
Lread_28_disablecopyonread_adam_m_dense_model_42_batch_normalization_85_gamma:Z
Lread_29_disablecopyonread_adam_v_dense_model_42_batch_normalization_85_gamma:Y
Kread_30_disablecopyonread_adam_m_dense_model_42_batch_normalization_85_beta:Y
Kread_31_disablecopyonread_adam_v_dense_model_42_batch_normalization_85_beta:R
@read_32_disablecopyonread_adam_m_dense_model_42_dense_128_kernel:R
@read_33_disablecopyonread_adam_v_dense_model_42_dense_128_kernel:L
>read_34_disablecopyonread_adam_m_dense_model_42_dense_128_bias:L
>read_35_disablecopyonread_adam_v_dense_model_42_dense_128_bias:+
!read_36_disablecopyonread_total_1: +
!read_37_disablecopyonread_count_1: )
read_38_disablecopyonread_total: )
read_39_disablecopyonread_count: 
savev2_const
identity_81��MergeV2Checkpoints�Read/DisableCopyOnRead�Read/ReadVariableOp�Read_1/DisableCopyOnRead�Read_1/ReadVariableOp�Read_10/DisableCopyOnRead�Read_10/ReadVariableOp�Read_11/DisableCopyOnRead�Read_11/ReadVariableOp�Read_12/DisableCopyOnRead�Read_12/ReadVariableOp�Read_13/DisableCopyOnRead�Read_13/ReadVariableOp�Read_14/DisableCopyOnRead�Read_14/ReadVariableOp�Read_15/DisableCopyOnRead�Read_15/ReadVariableOp�Read_16/DisableCopyOnRead�Read_16/ReadVariableOp�Read_17/DisableCopyOnRead�Read_17/ReadVariableOp�Read_18/DisableCopyOnRead�Read_18/ReadVariableOp�Read_19/DisableCopyOnRead�Read_19/ReadVariableOp�Read_2/DisableCopyOnRead�Read_2/ReadVariableOp�Read_20/DisableCopyOnRead�Read_20/ReadVariableOp�Read_21/DisableCopyOnRead�Read_21/ReadVariableOp�Read_22/DisableCopyOnRead�Read_22/ReadVariableOp�Read_23/DisableCopyOnRead�Read_23/ReadVariableOp�Read_24/DisableCopyOnRead�Read_24/ReadVariableOp�Read_25/DisableCopyOnRead�Read_25/ReadVariableOp�Read_26/DisableCopyOnRead�Read_26/ReadVariableOp�Read_27/DisableCopyOnRead�Read_27/ReadVariableOp�Read_28/DisableCopyOnRead�Read_28/ReadVariableOp�Read_29/DisableCopyOnRead�Read_29/ReadVariableOp�Read_3/DisableCopyOnRead�Read_3/ReadVariableOp�Read_30/DisableCopyOnRead�Read_30/ReadVariableOp�Read_31/DisableCopyOnRead�Read_31/ReadVariableOp�Read_32/DisableCopyOnRead�Read_32/ReadVariableOp�Read_33/DisableCopyOnRead�Read_33/ReadVariableOp�Read_34/DisableCopyOnRead�Read_34/ReadVariableOp�Read_35/DisableCopyOnRead�Read_35/ReadVariableOp�Read_36/DisableCopyOnRead�Read_36/ReadVariableOp�Read_37/DisableCopyOnRead�Read_37/ReadVariableOp�Read_38/DisableCopyOnRead�Read_38/ReadVariableOp�Read_39/DisableCopyOnRead�Read_39/ReadVariableOp�Read_4/DisableCopyOnRead�Read_4/ReadVariableOp�Read_5/DisableCopyOnRead�Read_5/ReadVariableOp�Read_6/DisableCopyOnRead�Read_6/ReadVariableOp�Read_7/DisableCopyOnRead�Read_7/ReadVariableOp�Read_8/DisableCopyOnRead�Read_8/ReadVariableOp�Read_9/DisableCopyOnRead�Read_9/ReadVariableOpw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
Read/DisableCopyOnReadDisableCopyOnRead6read_disablecopyonread_dense_model_42_dense_126_kernel"/device:CPU:0*
_output_shapes
 �
Read/ReadVariableOpReadVariableOp6read_disablecopyonread_dense_model_42_dense_126_kernel^Read/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0j
IdentityIdentityRead/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�b

Identity_1IdentityIdentity:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_1/DisableCopyOnReadDisableCopyOnRead6read_1_disablecopyonread_dense_model_42_dense_126_bias"/device:CPU:0*
_output_shapes
 �
Read_1/ReadVariableOpReadVariableOp6read_1_disablecopyonread_dense_model_42_dense_126_bias^Read_1/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_2IdentityRead_1/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_3IdentityIdentity_2:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_2/DisableCopyOnReadDisableCopyOnReadDread_2_disablecopyonread_dense_model_42_batch_normalization_84_gamma"/device:CPU:0*
_output_shapes
 �
Read_2/ReadVariableOpReadVariableOpDread_2_disablecopyonread_dense_model_42_batch_normalization_84_gamma^Read_2/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_4IdentityRead_2/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_5IdentityIdentity_4:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_3/DisableCopyOnReadDisableCopyOnReadCread_3_disablecopyonread_dense_model_42_batch_normalization_84_beta"/device:CPU:0*
_output_shapes
 �
Read_3/ReadVariableOpReadVariableOpCread_3_disablecopyonread_dense_model_42_batch_normalization_84_beta^Read_3/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_6IdentityRead_3/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_7IdentityIdentity_6:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_4/DisableCopyOnReadDisableCopyOnReadJread_4_disablecopyonread_dense_model_42_batch_normalization_84_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_4/ReadVariableOpReadVariableOpJread_4_disablecopyonread_dense_model_42_batch_normalization_84_moving_mean^Read_4/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_8IdentityRead_4/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_9IdentityIdentity_8:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_5/DisableCopyOnReadDisableCopyOnReadNread_5_disablecopyonread_dense_model_42_batch_normalization_84_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_5/ReadVariableOpReadVariableOpNread_5_disablecopyonread_dense_model_42_batch_normalization_84_moving_variance^Read_5/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_10IdentityRead_5/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_11IdentityIdentity_10:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_6/DisableCopyOnReadDisableCopyOnRead8read_6_disablecopyonread_dense_model_42_dense_127_kernel"/device:CPU:0*
_output_shapes
 �
Read_6/ReadVariableOpReadVariableOp8read_6_disablecopyonread_dense_model_42_dense_127_kernel^Read_6/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0n
Identity_12IdentityRead_6/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_13IdentityIdentity_12:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_7/DisableCopyOnReadDisableCopyOnRead6read_7_disablecopyonread_dense_model_42_dense_127_bias"/device:CPU:0*
_output_shapes
 �
Read_7/ReadVariableOpReadVariableOp6read_7_disablecopyonread_dense_model_42_dense_127_bias^Read_7/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_14IdentityRead_7/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_15IdentityIdentity_14:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_8/DisableCopyOnReadDisableCopyOnReadDread_8_disablecopyonread_dense_model_42_batch_normalization_85_gamma"/device:CPU:0*
_output_shapes
 �
Read_8/ReadVariableOpReadVariableOpDread_8_disablecopyonread_dense_model_42_batch_normalization_85_gamma^Read_8/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_16IdentityRead_8/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_17IdentityIdentity_16:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_9/DisableCopyOnReadDisableCopyOnReadCread_9_disablecopyonread_dense_model_42_batch_normalization_85_beta"/device:CPU:0*
_output_shapes
 �
Read_9/ReadVariableOpReadVariableOpCread_9_disablecopyonread_dense_model_42_batch_normalization_85_beta^Read_9/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_18IdentityRead_9/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_19IdentityIdentity_18:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_10/DisableCopyOnReadDisableCopyOnReadKread_10_disablecopyonread_dense_model_42_batch_normalization_85_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_10/ReadVariableOpReadVariableOpKread_10_disablecopyonread_dense_model_42_batch_normalization_85_moving_mean^Read_10/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_20IdentityRead_10/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_21IdentityIdentity_20:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_11/DisableCopyOnReadDisableCopyOnReadOread_11_disablecopyonread_dense_model_42_batch_normalization_85_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_11/ReadVariableOpReadVariableOpOread_11_disablecopyonread_dense_model_42_batch_normalization_85_moving_variance^Read_11/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_22IdentityRead_11/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_23IdentityIdentity_22:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_12/DisableCopyOnReadDisableCopyOnRead9read_12_disablecopyonread_dense_model_42_dense_128_kernel"/device:CPU:0*
_output_shapes
 �
Read_12/ReadVariableOpReadVariableOp9read_12_disablecopyonread_dense_model_42_dense_128_kernel^Read_12/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_24IdentityRead_12/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_25IdentityIdentity_24:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_13/DisableCopyOnReadDisableCopyOnRead7read_13_disablecopyonread_dense_model_42_dense_128_bias"/device:CPU:0*
_output_shapes
 �
Read_13/ReadVariableOpReadVariableOp7read_13_disablecopyonread_dense_model_42_dense_128_bias^Read_13/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_26IdentityRead_13/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_27IdentityIdentity_26:output:0"/device:CPU:0*
T0*
_output_shapes
:x
Read_14/DisableCopyOnReadDisableCopyOnRead#read_14_disablecopyonread_iteration"/device:CPU:0*
_output_shapes
 �
Read_14/ReadVariableOpReadVariableOp#read_14_disablecopyonread_iteration^Read_14/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0	g
Identity_28IdentityRead_14/ReadVariableOp:value:0"/device:CPU:0*
T0	*
_output_shapes
: ]
Identity_29IdentityIdentity_28:output:0"/device:CPU:0*
T0	*
_output_shapes
: |
Read_15/DisableCopyOnReadDisableCopyOnRead'read_15_disablecopyonread_learning_rate"/device:CPU:0*
_output_shapes
 �
Read_15/ReadVariableOpReadVariableOp'read_15_disablecopyonread_learning_rate^Read_15/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_30IdentityRead_15/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_31IdentityIdentity_30:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_16/DisableCopyOnReadDisableCopyOnRead@read_16_disablecopyonread_adam_m_dense_model_42_dense_126_kernel"/device:CPU:0*
_output_shapes
 �
Read_16/ReadVariableOpReadVariableOp@read_16_disablecopyonread_adam_m_dense_model_42_dense_126_kernel^Read_16/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0p
Identity_32IdentityRead_16/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�f
Identity_33IdentityIdentity_32:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_17/DisableCopyOnReadDisableCopyOnRead@read_17_disablecopyonread_adam_v_dense_model_42_dense_126_kernel"/device:CPU:0*
_output_shapes
 �
Read_17/ReadVariableOpReadVariableOp@read_17_disablecopyonread_adam_v_dense_model_42_dense_126_kernel^Read_17/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0p
Identity_34IdentityRead_17/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�f
Identity_35IdentityIdentity_34:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_18/DisableCopyOnReadDisableCopyOnRead>read_18_disablecopyonread_adam_m_dense_model_42_dense_126_bias"/device:CPU:0*
_output_shapes
 �
Read_18/ReadVariableOpReadVariableOp>read_18_disablecopyonread_adam_m_dense_model_42_dense_126_bias^Read_18/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_36IdentityRead_18/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_37IdentityIdentity_36:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_19/DisableCopyOnReadDisableCopyOnRead>read_19_disablecopyonread_adam_v_dense_model_42_dense_126_bias"/device:CPU:0*
_output_shapes
 �
Read_19/ReadVariableOpReadVariableOp>read_19_disablecopyonread_adam_v_dense_model_42_dense_126_bias^Read_19/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_38IdentityRead_19/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_39IdentityIdentity_38:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_20/DisableCopyOnReadDisableCopyOnReadLread_20_disablecopyonread_adam_m_dense_model_42_batch_normalization_84_gamma"/device:CPU:0*
_output_shapes
 �
Read_20/ReadVariableOpReadVariableOpLread_20_disablecopyonread_adam_m_dense_model_42_batch_normalization_84_gamma^Read_20/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_40IdentityRead_20/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_41IdentityIdentity_40:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_21/DisableCopyOnReadDisableCopyOnReadLread_21_disablecopyonread_adam_v_dense_model_42_batch_normalization_84_gamma"/device:CPU:0*
_output_shapes
 �
Read_21/ReadVariableOpReadVariableOpLread_21_disablecopyonread_adam_v_dense_model_42_batch_normalization_84_gamma^Read_21/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_42IdentityRead_21/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_43IdentityIdentity_42:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_22/DisableCopyOnReadDisableCopyOnReadKread_22_disablecopyonread_adam_m_dense_model_42_batch_normalization_84_beta"/device:CPU:0*
_output_shapes
 �
Read_22/ReadVariableOpReadVariableOpKread_22_disablecopyonread_adam_m_dense_model_42_batch_normalization_84_beta^Read_22/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_44IdentityRead_22/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_45IdentityIdentity_44:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_23/DisableCopyOnReadDisableCopyOnReadKread_23_disablecopyonread_adam_v_dense_model_42_batch_normalization_84_beta"/device:CPU:0*
_output_shapes
 �
Read_23/ReadVariableOpReadVariableOpKread_23_disablecopyonread_adam_v_dense_model_42_batch_normalization_84_beta^Read_23/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_46IdentityRead_23/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_47IdentityIdentity_46:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_24/DisableCopyOnReadDisableCopyOnRead@read_24_disablecopyonread_adam_m_dense_model_42_dense_127_kernel"/device:CPU:0*
_output_shapes
 �
Read_24/ReadVariableOpReadVariableOp@read_24_disablecopyonread_adam_m_dense_model_42_dense_127_kernel^Read_24/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_48IdentityRead_24/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_49IdentityIdentity_48:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_25/DisableCopyOnReadDisableCopyOnRead@read_25_disablecopyonread_adam_v_dense_model_42_dense_127_kernel"/device:CPU:0*
_output_shapes
 �
Read_25/ReadVariableOpReadVariableOp@read_25_disablecopyonread_adam_v_dense_model_42_dense_127_kernel^Read_25/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_50IdentityRead_25/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_51IdentityIdentity_50:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_26/DisableCopyOnReadDisableCopyOnRead>read_26_disablecopyonread_adam_m_dense_model_42_dense_127_bias"/device:CPU:0*
_output_shapes
 �
Read_26/ReadVariableOpReadVariableOp>read_26_disablecopyonread_adam_m_dense_model_42_dense_127_bias^Read_26/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_52IdentityRead_26/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_53IdentityIdentity_52:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_27/DisableCopyOnReadDisableCopyOnRead>read_27_disablecopyonread_adam_v_dense_model_42_dense_127_bias"/device:CPU:0*
_output_shapes
 �
Read_27/ReadVariableOpReadVariableOp>read_27_disablecopyonread_adam_v_dense_model_42_dense_127_bias^Read_27/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_54IdentityRead_27/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_55IdentityIdentity_54:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_28/DisableCopyOnReadDisableCopyOnReadLread_28_disablecopyonread_adam_m_dense_model_42_batch_normalization_85_gamma"/device:CPU:0*
_output_shapes
 �
Read_28/ReadVariableOpReadVariableOpLread_28_disablecopyonread_adam_m_dense_model_42_batch_normalization_85_gamma^Read_28/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_56IdentityRead_28/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_57IdentityIdentity_56:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_29/DisableCopyOnReadDisableCopyOnReadLread_29_disablecopyonread_adam_v_dense_model_42_batch_normalization_85_gamma"/device:CPU:0*
_output_shapes
 �
Read_29/ReadVariableOpReadVariableOpLread_29_disablecopyonread_adam_v_dense_model_42_batch_normalization_85_gamma^Read_29/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_58IdentityRead_29/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_59IdentityIdentity_58:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_30/DisableCopyOnReadDisableCopyOnReadKread_30_disablecopyonread_adam_m_dense_model_42_batch_normalization_85_beta"/device:CPU:0*
_output_shapes
 �
Read_30/ReadVariableOpReadVariableOpKread_30_disablecopyonread_adam_m_dense_model_42_batch_normalization_85_beta^Read_30/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_60IdentityRead_30/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_61IdentityIdentity_60:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_31/DisableCopyOnReadDisableCopyOnReadKread_31_disablecopyonread_adam_v_dense_model_42_batch_normalization_85_beta"/device:CPU:0*
_output_shapes
 �
Read_31/ReadVariableOpReadVariableOpKread_31_disablecopyonread_adam_v_dense_model_42_batch_normalization_85_beta^Read_31/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_62IdentityRead_31/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_63IdentityIdentity_62:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_32/DisableCopyOnReadDisableCopyOnRead@read_32_disablecopyonread_adam_m_dense_model_42_dense_128_kernel"/device:CPU:0*
_output_shapes
 �
Read_32/ReadVariableOpReadVariableOp@read_32_disablecopyonread_adam_m_dense_model_42_dense_128_kernel^Read_32/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_64IdentityRead_32/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_65IdentityIdentity_64:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_33/DisableCopyOnReadDisableCopyOnRead@read_33_disablecopyonread_adam_v_dense_model_42_dense_128_kernel"/device:CPU:0*
_output_shapes
 �
Read_33/ReadVariableOpReadVariableOp@read_33_disablecopyonread_adam_v_dense_model_42_dense_128_kernel^Read_33/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_66IdentityRead_33/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_67IdentityIdentity_66:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_34/DisableCopyOnReadDisableCopyOnRead>read_34_disablecopyonread_adam_m_dense_model_42_dense_128_bias"/device:CPU:0*
_output_shapes
 �
Read_34/ReadVariableOpReadVariableOp>read_34_disablecopyonread_adam_m_dense_model_42_dense_128_bias^Read_34/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_68IdentityRead_34/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_69IdentityIdentity_68:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_35/DisableCopyOnReadDisableCopyOnRead>read_35_disablecopyonread_adam_v_dense_model_42_dense_128_bias"/device:CPU:0*
_output_shapes
 �
Read_35/ReadVariableOpReadVariableOp>read_35_disablecopyonread_adam_v_dense_model_42_dense_128_bias^Read_35/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_70IdentityRead_35/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_71IdentityIdentity_70:output:0"/device:CPU:0*
T0*
_output_shapes
:v
Read_36/DisableCopyOnReadDisableCopyOnRead!read_36_disablecopyonread_total_1"/device:CPU:0*
_output_shapes
 �
Read_36/ReadVariableOpReadVariableOp!read_36_disablecopyonread_total_1^Read_36/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_72IdentityRead_36/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_73IdentityIdentity_72:output:0"/device:CPU:0*
T0*
_output_shapes
: v
Read_37/DisableCopyOnReadDisableCopyOnRead!read_37_disablecopyonread_count_1"/device:CPU:0*
_output_shapes
 �
Read_37/ReadVariableOpReadVariableOp!read_37_disablecopyonread_count_1^Read_37/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_74IdentityRead_37/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_75IdentityIdentity_74:output:0"/device:CPU:0*
T0*
_output_shapes
: t
Read_38/DisableCopyOnReadDisableCopyOnReadread_38_disablecopyonread_total"/device:CPU:0*
_output_shapes
 �
Read_38/ReadVariableOpReadVariableOpread_38_disablecopyonread_total^Read_38/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_76IdentityRead_38/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_77IdentityIdentity_76:output:0"/device:CPU:0*
T0*
_output_shapes
: t
Read_39/DisableCopyOnReadDisableCopyOnReadread_39_disablecopyonread_count"/device:CPU:0*
_output_shapes
 �
Read_39/ReadVariableOpReadVariableOpread_39_disablecopyonread_count^Read_39/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_78IdentityRead_39/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_79IdentityIdentity_78:output:0"/device:CPU:0*
T0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:)*
dtype0*�
value�B�)B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:)*
dtype0*e
value\BZ)B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0Identity_1:output:0Identity_3:output:0Identity_5:output:0Identity_7:output:0Identity_9:output:0Identity_11:output:0Identity_13:output:0Identity_15:output:0Identity_17:output:0Identity_19:output:0Identity_21:output:0Identity_23:output:0Identity_25:output:0Identity_27:output:0Identity_29:output:0Identity_31:output:0Identity_33:output:0Identity_35:output:0Identity_37:output:0Identity_39:output:0Identity_41:output:0Identity_43:output:0Identity_45:output:0Identity_47:output:0Identity_49:output:0Identity_51:output:0Identity_53:output:0Identity_55:output:0Identity_57:output:0Identity_59:output:0Identity_61:output:0Identity_63:output:0Identity_65:output:0Identity_67:output:0Identity_69:output:0Identity_71:output:0Identity_73:output:0Identity_75:output:0Identity_77:output:0Identity_79:output:0savev2_const"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *7
dtypes-
+2)	�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 i
Identity_80Identityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: U
Identity_81IdentityIdentity_80:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^MergeV2Checkpoints^Read/DisableCopyOnRead^Read/ReadVariableOp^Read_1/DisableCopyOnRead^Read_1/ReadVariableOp^Read_10/DisableCopyOnRead^Read_10/ReadVariableOp^Read_11/DisableCopyOnRead^Read_11/ReadVariableOp^Read_12/DisableCopyOnRead^Read_12/ReadVariableOp^Read_13/DisableCopyOnRead^Read_13/ReadVariableOp^Read_14/DisableCopyOnRead^Read_14/ReadVariableOp^Read_15/DisableCopyOnRead^Read_15/ReadVariableOp^Read_16/DisableCopyOnRead^Read_16/ReadVariableOp^Read_17/DisableCopyOnRead^Read_17/ReadVariableOp^Read_18/DisableCopyOnRead^Read_18/ReadVariableOp^Read_19/DisableCopyOnRead^Read_19/ReadVariableOp^Read_2/DisableCopyOnRead^Read_2/ReadVariableOp^Read_20/DisableCopyOnRead^Read_20/ReadVariableOp^Read_21/DisableCopyOnRead^Read_21/ReadVariableOp^Read_22/DisableCopyOnRead^Read_22/ReadVariableOp^Read_23/DisableCopyOnRead^Read_23/ReadVariableOp^Read_24/DisableCopyOnRead^Read_24/ReadVariableOp^Read_25/DisableCopyOnRead^Read_25/ReadVariableOp^Read_26/DisableCopyOnRead^Read_26/ReadVariableOp^Read_27/DisableCopyOnRead^Read_27/ReadVariableOp^Read_28/DisableCopyOnRead^Read_28/ReadVariableOp^Read_29/DisableCopyOnRead^Read_29/ReadVariableOp^Read_3/DisableCopyOnRead^Read_3/ReadVariableOp^Read_30/DisableCopyOnRead^Read_30/ReadVariableOp^Read_31/DisableCopyOnRead^Read_31/ReadVariableOp^Read_32/DisableCopyOnRead^Read_32/ReadVariableOp^Read_33/DisableCopyOnRead^Read_33/ReadVariableOp^Read_34/DisableCopyOnRead^Read_34/ReadVariableOp^Read_35/DisableCopyOnRead^Read_35/ReadVariableOp^Read_36/DisableCopyOnRead^Read_36/ReadVariableOp^Read_37/DisableCopyOnRead^Read_37/ReadVariableOp^Read_38/DisableCopyOnRead^Read_38/ReadVariableOp^Read_39/DisableCopyOnRead^Read_39/ReadVariableOp^Read_4/DisableCopyOnRead^Read_4/ReadVariableOp^Read_5/DisableCopyOnRead^Read_5/ReadVariableOp^Read_6/DisableCopyOnRead^Read_6/ReadVariableOp^Read_7/DisableCopyOnRead^Read_7/ReadVariableOp^Read_8/DisableCopyOnRead^Read_8/ReadVariableOp^Read_9/DisableCopyOnRead^Read_9/ReadVariableOp*
_output_shapes
 "#
identity_81Identity_81:output:0*(
_construction_contextkEagerRuntime*g
_input_shapesV
T: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2(
MergeV2CheckpointsMergeV2Checkpoints20
Read/DisableCopyOnReadRead/DisableCopyOnRead2*
Read/ReadVariableOpRead/ReadVariableOp24
Read_1/DisableCopyOnReadRead_1/DisableCopyOnRead2.
Read_1/ReadVariableOpRead_1/ReadVariableOp26
Read_10/DisableCopyOnReadRead_10/DisableCopyOnRead20
Read_10/ReadVariableOpRead_10/ReadVariableOp26
Read_11/DisableCopyOnReadRead_11/DisableCopyOnRead20
Read_11/ReadVariableOpRead_11/ReadVariableOp26
Read_12/DisableCopyOnReadRead_12/DisableCopyOnRead20
Read_12/ReadVariableOpRead_12/ReadVariableOp26
Read_13/DisableCopyOnReadRead_13/DisableCopyOnRead20
Read_13/ReadVariableOpRead_13/ReadVariableOp26
Read_14/DisableCopyOnReadRead_14/DisableCopyOnRead20
Read_14/ReadVariableOpRead_14/ReadVariableOp26
Read_15/DisableCopyOnReadRead_15/DisableCopyOnRead20
Read_15/ReadVariableOpRead_15/ReadVariableOp26
Read_16/DisableCopyOnReadRead_16/DisableCopyOnRead20
Read_16/ReadVariableOpRead_16/ReadVariableOp26
Read_17/DisableCopyOnReadRead_17/DisableCopyOnRead20
Read_17/ReadVariableOpRead_17/ReadVariableOp26
Read_18/DisableCopyOnReadRead_18/DisableCopyOnRead20
Read_18/ReadVariableOpRead_18/ReadVariableOp26
Read_19/DisableCopyOnReadRead_19/DisableCopyOnRead20
Read_19/ReadVariableOpRead_19/ReadVariableOp24
Read_2/DisableCopyOnReadRead_2/DisableCopyOnRead2.
Read_2/ReadVariableOpRead_2/ReadVariableOp26
Read_20/DisableCopyOnReadRead_20/DisableCopyOnRead20
Read_20/ReadVariableOpRead_20/ReadVariableOp26
Read_21/DisableCopyOnReadRead_21/DisableCopyOnRead20
Read_21/ReadVariableOpRead_21/ReadVariableOp26
Read_22/DisableCopyOnReadRead_22/DisableCopyOnRead20
Read_22/ReadVariableOpRead_22/ReadVariableOp26
Read_23/DisableCopyOnReadRead_23/DisableCopyOnRead20
Read_23/ReadVariableOpRead_23/ReadVariableOp26
Read_24/DisableCopyOnReadRead_24/DisableCopyOnRead20
Read_24/ReadVariableOpRead_24/ReadVariableOp26
Read_25/DisableCopyOnReadRead_25/DisableCopyOnRead20
Read_25/ReadVariableOpRead_25/ReadVariableOp26
Read_26/DisableCopyOnReadRead_26/DisableCopyOnRead20
Read_26/ReadVariableOpRead_26/ReadVariableOp26
Read_27/DisableCopyOnReadRead_27/DisableCopyOnRead20
Read_27/ReadVariableOpRead_27/ReadVariableOp26
Read_28/DisableCopyOnReadRead_28/DisableCopyOnRead20
Read_28/ReadVariableOpRead_28/ReadVariableOp26
Read_29/DisableCopyOnReadRead_29/DisableCopyOnRead20
Read_29/ReadVariableOpRead_29/ReadVariableOp24
Read_3/DisableCopyOnReadRead_3/DisableCopyOnRead2.
Read_3/ReadVariableOpRead_3/ReadVariableOp26
Read_30/DisableCopyOnReadRead_30/DisableCopyOnRead20
Read_30/ReadVariableOpRead_30/ReadVariableOp26
Read_31/DisableCopyOnReadRead_31/DisableCopyOnRead20
Read_31/ReadVariableOpRead_31/ReadVariableOp26
Read_32/DisableCopyOnReadRead_32/DisableCopyOnRead20
Read_32/ReadVariableOpRead_32/ReadVariableOp26
Read_33/DisableCopyOnReadRead_33/DisableCopyOnRead20
Read_33/ReadVariableOpRead_33/ReadVariableOp26
Read_34/DisableCopyOnReadRead_34/DisableCopyOnRead20
Read_34/ReadVariableOpRead_34/ReadVariableOp26
Read_35/DisableCopyOnReadRead_35/DisableCopyOnRead20
Read_35/ReadVariableOpRead_35/ReadVariableOp26
Read_36/DisableCopyOnReadRead_36/DisableCopyOnRead20
Read_36/ReadVariableOpRead_36/ReadVariableOp26
Read_37/DisableCopyOnReadRead_37/DisableCopyOnRead20
Read_37/ReadVariableOpRead_37/ReadVariableOp26
Read_38/DisableCopyOnReadRead_38/DisableCopyOnRead20
Read_38/ReadVariableOpRead_38/ReadVariableOp26
Read_39/DisableCopyOnReadRead_39/DisableCopyOnRead20
Read_39/ReadVariableOpRead_39/ReadVariableOp24
Read_4/DisableCopyOnReadRead_4/DisableCopyOnRead2.
Read_4/ReadVariableOpRead_4/ReadVariableOp24
Read_5/DisableCopyOnReadRead_5/DisableCopyOnRead2.
Read_5/ReadVariableOpRead_5/ReadVariableOp24
Read_6/DisableCopyOnReadRead_6/DisableCopyOnRead2.
Read_6/ReadVariableOpRead_6/ReadVariableOp24
Read_7/DisableCopyOnReadRead_7/DisableCopyOnRead2.
Read_7/ReadVariableOpRead_7/ReadVariableOp24
Read_8/DisableCopyOnReadRead_8/DisableCopyOnRead2.
Read_8/ReadVariableOpRead_8/ReadVariableOp24
Read_9/DisableCopyOnReadRead_9/DisableCopyOnRead2.
Read_9/ReadVariableOpRead_9/ReadVariableOp:=)9

_output_shapes
: 

_user_specified_nameConst:%(!

_user_specified_namecount:%'!

_user_specified_nametotal:'&#
!
_user_specified_name	count_1:'%#
!
_user_specified_name	total_1:D$@
>
_user_specified_name&$Adam/v/dense_model_42/dense_128/bias:D#@
>
_user_specified_name&$Adam/m/dense_model_42/dense_128/bias:F"B
@
_user_specified_name(&Adam/v/dense_model_42/dense_128/kernel:F!B
@
_user_specified_name(&Adam/m/dense_model_42/dense_128/kernel:Q M
K
_user_specified_name31Adam/v/dense_model_42/batch_normalization_85/beta:QM
K
_user_specified_name31Adam/m/dense_model_42/batch_normalization_85/beta:RN
L
_user_specified_name42Adam/v/dense_model_42/batch_normalization_85/gamma:RN
L
_user_specified_name42Adam/m/dense_model_42/batch_normalization_85/gamma:D@
>
_user_specified_name&$Adam/v/dense_model_42/dense_127/bias:D@
>
_user_specified_name&$Adam/m/dense_model_42/dense_127/bias:FB
@
_user_specified_name(&Adam/v/dense_model_42/dense_127/kernel:FB
@
_user_specified_name(&Adam/m/dense_model_42/dense_127/kernel:QM
K
_user_specified_name31Adam/v/dense_model_42/batch_normalization_84/beta:QM
K
_user_specified_name31Adam/m/dense_model_42/batch_normalization_84/beta:RN
L
_user_specified_name42Adam/v/dense_model_42/batch_normalization_84/gamma:RN
L
_user_specified_name42Adam/m/dense_model_42/batch_normalization_84/gamma:D@
>
_user_specified_name&$Adam/v/dense_model_42/dense_126/bias:D@
>
_user_specified_name&$Adam/m/dense_model_42/dense_126/bias:FB
@
_user_specified_name(&Adam/v/dense_model_42/dense_126/kernel:FB
@
_user_specified_name(&Adam/m/dense_model_42/dense_126/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:=9
7
_user_specified_namedense_model_42/dense_128/bias:?;
9
_user_specified_name!dense_model_42/dense_128/kernel:UQ
O
_user_specified_name75dense_model_42/batch_normalization_85/moving_variance:QM
K
_user_specified_name31dense_model_42/batch_normalization_85/moving_mean:J
F
D
_user_specified_name,*dense_model_42/batch_normalization_85/beta:K	G
E
_user_specified_name-+dense_model_42/batch_normalization_85/gamma:=9
7
_user_specified_namedense_model_42/dense_127/bias:?;
9
_user_specified_name!dense_model_42/dense_127/kernel:UQ
O
_user_specified_name75dense_model_42/batch_normalization_84/moving_variance:QM
K
_user_specified_name31dense_model_42/batch_normalization_84/moving_mean:JF
D
_user_specified_name,*dense_model_42/batch_normalization_84/beta:KG
E
_user_specified_name-+dense_model_42/batch_normalization_84/gamma:=9
7
_user_specified_namedense_model_42/dense_126/bias:?;
9
_user_specified_name!dense_model_42/dense_126/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
$__inference_signature_wrapper_161720
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference__wrapped_model_161326o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161716:&"
 
_user_specified_name161714:&"
 
_user_specified_name161712:&"
 
_user_specified_name161710:&
"
 
_user_specified_name161708:&	"
 
_user_specified_name161706:&"
 
_user_specified_name161704:&"
 
_user_specified_name161702:&"
 
_user_specified_name161700:&"
 
_user_specified_name161698:&"
 
_user_specified_name161696:&"
 
_user_specified_name161694:&"
 
_user_specified_name161692:&"
 
_user_specified_name161690:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�`
�
!__inference__wrapped_model_161326
input_1J
7dense_model_42_dense_126_matmul_readvariableop_resource:	�F
8dense_model_42_dense_126_biasadd_readvariableop_resource:P
Bdense_model_42_batch_normalization_84_cast_readvariableop_resource:R
Ddense_model_42_batch_normalization_84_cast_1_readvariableop_resource:R
Ddense_model_42_batch_normalization_84_cast_2_readvariableop_resource:R
Ddense_model_42_batch_normalization_84_cast_3_readvariableop_resource:I
7dense_model_42_dense_127_matmul_readvariableop_resource:F
8dense_model_42_dense_127_biasadd_readvariableop_resource:P
Bdense_model_42_batch_normalization_85_cast_readvariableop_resource:R
Ddense_model_42_batch_normalization_85_cast_1_readvariableop_resource:R
Ddense_model_42_batch_normalization_85_cast_2_readvariableop_resource:R
Ddense_model_42_batch_normalization_85_cast_3_readvariableop_resource:I
7dense_model_42_dense_128_matmul_readvariableop_resource:F
8dense_model_42_dense_128_biasadd_readvariableop_resource:
identity��9dense_model_42/batch_normalization_84/Cast/ReadVariableOp�;dense_model_42/batch_normalization_84/Cast_1/ReadVariableOp�;dense_model_42/batch_normalization_84/Cast_2/ReadVariableOp�;dense_model_42/batch_normalization_84/Cast_3/ReadVariableOp�9dense_model_42/batch_normalization_85/Cast/ReadVariableOp�;dense_model_42/batch_normalization_85/Cast_1/ReadVariableOp�;dense_model_42/batch_normalization_85/Cast_2/ReadVariableOp�;dense_model_42/batch_normalization_85/Cast_3/ReadVariableOp�/dense_model_42/dense_126/BiasAdd/ReadVariableOp�.dense_model_42/dense_126/MatMul/ReadVariableOp�/dense_model_42/dense_127/BiasAdd/ReadVariableOp�.dense_model_42/dense_127/MatMul/ReadVariableOp�/dense_model_42/dense_128/BiasAdd/ReadVariableOp�.dense_model_42/dense_128/MatMul/ReadVariableOp�
.dense_model_42/dense_126/MatMul/ReadVariableOpReadVariableOp7dense_model_42_dense_126_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_model_42/dense_126/MatMulMatMulinput_16dense_model_42/dense_126/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
/dense_model_42/dense_126/BiasAdd/ReadVariableOpReadVariableOp8dense_model_42_dense_126_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_model_42/dense_126/BiasAddBiasAdd)dense_model_42/dense_126/MatMul:product:07dense_model_42/dense_126/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
9dense_model_42/batch_normalization_84/Cast/ReadVariableOpReadVariableOpBdense_model_42_batch_normalization_84_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
;dense_model_42/batch_normalization_84/Cast_1/ReadVariableOpReadVariableOpDdense_model_42_batch_normalization_84_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
;dense_model_42/batch_normalization_84/Cast_2/ReadVariableOpReadVariableOpDdense_model_42_batch_normalization_84_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
;dense_model_42/batch_normalization_84/Cast_3/ReadVariableOpReadVariableOpDdense_model_42_batch_normalization_84_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0z
5dense_model_42/batch_normalization_84/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
3dense_model_42/batch_normalization_84/batchnorm/addAddV2Cdense_model_42/batch_normalization_84/Cast_1/ReadVariableOp:value:0>dense_model_42/batch_normalization_84/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
5dense_model_42/batch_normalization_84/batchnorm/RsqrtRsqrt7dense_model_42/batch_normalization_84/batchnorm/add:z:0*
T0*
_output_shapes
:�
3dense_model_42/batch_normalization_84/batchnorm/mulMul9dense_model_42/batch_normalization_84/batchnorm/Rsqrt:y:0Cdense_model_42/batch_normalization_84/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
5dense_model_42/batch_normalization_84/batchnorm/mul_1Mul)dense_model_42/dense_126/BiasAdd:output:07dense_model_42/batch_normalization_84/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
5dense_model_42/batch_normalization_84/batchnorm/mul_2MulAdense_model_42/batch_normalization_84/Cast/ReadVariableOp:value:07dense_model_42/batch_normalization_84/batchnorm/mul:z:0*
T0*
_output_shapes
:�
3dense_model_42/batch_normalization_84/batchnorm/subSubCdense_model_42/batch_normalization_84/Cast_2/ReadVariableOp:value:09dense_model_42/batch_normalization_84/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
5dense_model_42/batch_normalization_84/batchnorm/add_1AddV29dense_model_42/batch_normalization_84/batchnorm/mul_1:z:07dense_model_42/batch_normalization_84/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
.dense_model_42/dense_127/MatMul/ReadVariableOpReadVariableOp7dense_model_42_dense_127_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_model_42/dense_127/MatMulMatMul9dense_model_42/batch_normalization_84/batchnorm/add_1:z:06dense_model_42/dense_127/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
/dense_model_42/dense_127/BiasAdd/ReadVariableOpReadVariableOp8dense_model_42_dense_127_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_model_42/dense_127/BiasAddBiasAdd)dense_model_42/dense_127/MatMul:product:07dense_model_42/dense_127/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
9dense_model_42/batch_normalization_85/Cast/ReadVariableOpReadVariableOpBdense_model_42_batch_normalization_85_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
;dense_model_42/batch_normalization_85/Cast_1/ReadVariableOpReadVariableOpDdense_model_42_batch_normalization_85_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
;dense_model_42/batch_normalization_85/Cast_2/ReadVariableOpReadVariableOpDdense_model_42_batch_normalization_85_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
;dense_model_42/batch_normalization_85/Cast_3/ReadVariableOpReadVariableOpDdense_model_42_batch_normalization_85_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0z
5dense_model_42/batch_normalization_85/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
3dense_model_42/batch_normalization_85/batchnorm/addAddV2Cdense_model_42/batch_normalization_85/Cast_1/ReadVariableOp:value:0>dense_model_42/batch_normalization_85/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
5dense_model_42/batch_normalization_85/batchnorm/RsqrtRsqrt7dense_model_42/batch_normalization_85/batchnorm/add:z:0*
T0*
_output_shapes
:�
3dense_model_42/batch_normalization_85/batchnorm/mulMul9dense_model_42/batch_normalization_85/batchnorm/Rsqrt:y:0Cdense_model_42/batch_normalization_85/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
5dense_model_42/batch_normalization_85/batchnorm/mul_1Mul)dense_model_42/dense_127/BiasAdd:output:07dense_model_42/batch_normalization_85/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
5dense_model_42/batch_normalization_85/batchnorm/mul_2MulAdense_model_42/batch_normalization_85/Cast/ReadVariableOp:value:07dense_model_42/batch_normalization_85/batchnorm/mul:z:0*
T0*
_output_shapes
:�
3dense_model_42/batch_normalization_85/batchnorm/subSubCdense_model_42/batch_normalization_85/Cast_2/ReadVariableOp:value:09dense_model_42/batch_normalization_85/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
5dense_model_42/batch_normalization_85/batchnorm/add_1AddV29dense_model_42/batch_normalization_85/batchnorm/mul_1:z:07dense_model_42/batch_normalization_85/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
.dense_model_42/dense_128/MatMul/ReadVariableOpReadVariableOp7dense_model_42_dense_128_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_model_42/dense_128/MatMulMatMul9dense_model_42/batch_normalization_85/batchnorm/add_1:z:06dense_model_42/dense_128/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
/dense_model_42/dense_128/BiasAdd/ReadVariableOpReadVariableOp8dense_model_42_dense_128_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_model_42/dense_128/BiasAddBiasAdd)dense_model_42/dense_128/MatMul:product:07dense_model_42/dense_128/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 dense_model_42/dense_128/SigmoidSigmoid)dense_model_42/dense_128/BiasAdd:output:0*
T0*'
_output_shapes
:���������s
IdentityIdentity$dense_model_42/dense_128/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp:^dense_model_42/batch_normalization_84/Cast/ReadVariableOp<^dense_model_42/batch_normalization_84/Cast_1/ReadVariableOp<^dense_model_42/batch_normalization_84/Cast_2/ReadVariableOp<^dense_model_42/batch_normalization_84/Cast_3/ReadVariableOp:^dense_model_42/batch_normalization_85/Cast/ReadVariableOp<^dense_model_42/batch_normalization_85/Cast_1/ReadVariableOp<^dense_model_42/batch_normalization_85/Cast_2/ReadVariableOp<^dense_model_42/batch_normalization_85/Cast_3/ReadVariableOp0^dense_model_42/dense_126/BiasAdd/ReadVariableOp/^dense_model_42/dense_126/MatMul/ReadVariableOp0^dense_model_42/dense_127/BiasAdd/ReadVariableOp/^dense_model_42/dense_127/MatMul/ReadVariableOp0^dense_model_42/dense_128/BiasAdd/ReadVariableOp/^dense_model_42/dense_128/MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2v
9dense_model_42/batch_normalization_84/Cast/ReadVariableOp9dense_model_42/batch_normalization_84/Cast/ReadVariableOp2z
;dense_model_42/batch_normalization_84/Cast_1/ReadVariableOp;dense_model_42/batch_normalization_84/Cast_1/ReadVariableOp2z
;dense_model_42/batch_normalization_84/Cast_2/ReadVariableOp;dense_model_42/batch_normalization_84/Cast_2/ReadVariableOp2z
;dense_model_42/batch_normalization_84/Cast_3/ReadVariableOp;dense_model_42/batch_normalization_84/Cast_3/ReadVariableOp2v
9dense_model_42/batch_normalization_85/Cast/ReadVariableOp9dense_model_42/batch_normalization_85/Cast/ReadVariableOp2z
;dense_model_42/batch_normalization_85/Cast_1/ReadVariableOp;dense_model_42/batch_normalization_85/Cast_1/ReadVariableOp2z
;dense_model_42/batch_normalization_85/Cast_2/ReadVariableOp;dense_model_42/batch_normalization_85/Cast_2/ReadVariableOp2z
;dense_model_42/batch_normalization_85/Cast_3/ReadVariableOp;dense_model_42/batch_normalization_85/Cast_3/ReadVariableOp2b
/dense_model_42/dense_126/BiasAdd/ReadVariableOp/dense_model_42/dense_126/BiasAdd/ReadVariableOp2`
.dense_model_42/dense_126/MatMul/ReadVariableOp.dense_model_42/dense_126/MatMul/ReadVariableOp2b
/dense_model_42/dense_127/BiasAdd/ReadVariableOp/dense_model_42/dense_127/BiasAdd/ReadVariableOp2`
.dense_model_42/dense_127/MatMul/ReadVariableOp.dense_model_42/dense_127/MatMul/ReadVariableOp2b
/dense_model_42/dense_128/BiasAdd/ReadVariableOp/dense_model_42/dense_128/BiasAdd/ReadVariableOp2`
.dense_model_42/dense_128/MatMul/ReadVariableOp.dense_model_42/dense_128/MatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:(
$
"
_user_specified_name
resource:(	$
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�%
�
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161898

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 Z
AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
AssignMovingAvg/ReadVariableOpReadVariableOp'assignmovingavg_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvgAssignSubVariableOp'assignmovingavg_readvariableop_resourceAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0\
AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
 AssignMovingAvg_1/ReadVariableOpReadVariableOp)assignmovingavg_1_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:q
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12"
AssignMovingAvgAssignMovingAvg2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
E__inference_dense_128_layer_call_and_return_conditional_losses_161547

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:���������Z
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�$
�
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161554
input_1#
dense_126_161499:	�
dense_126_161501:+
batch_normalization_84_161504:+
batch_normalization_84_161506:+
batch_normalization_84_161508:+
batch_normalization_84_161510:"
dense_127_161523:
dense_127_161525:+
batch_normalization_85_161528:+
batch_normalization_85_161530:+
batch_normalization_85_161532:+
batch_normalization_85_161534:"
dense_128_161548:
dense_128_161550:
identity��.batch_normalization_84/StatefulPartitionedCall�.batch_normalization_85/StatefulPartitionedCall�!dense_126/StatefulPartitionedCall�!dense_127/StatefulPartitionedCall�!dense_128/StatefulPartitionedCall�
!dense_126/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_126_161499dense_126_161501*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_126_layer_call_and_return_conditional_losses_161498�
.batch_normalization_84/StatefulPartitionedCallStatefulPartitionedCall*dense_126/StatefulPartitionedCall:output:0batch_normalization_84_161504batch_normalization_84_161506batch_normalization_84_161508batch_normalization_84_161510*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161360�
!dense_127/StatefulPartitionedCallStatefulPartitionedCall7batch_normalization_84/StatefulPartitionedCall:output:0dense_127_161523dense_127_161525*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_127_layer_call_and_return_conditional_losses_161522�
.batch_normalization_85/StatefulPartitionedCallStatefulPartitionedCall*dense_127/StatefulPartitionedCall:output:0batch_normalization_85_161528batch_normalization_85_161530batch_normalization_85_161532batch_normalization_85_161534*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161440�
!dense_128/StatefulPartitionedCallStatefulPartitionedCall7batch_normalization_85/StatefulPartitionedCall:output:0dense_128_161548dense_128_161550*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_128_layer_call_and_return_conditional_losses_161547y
IdentityIdentity*dense_128/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp/^batch_normalization_84/StatefulPartitionedCall/^batch_normalization_85/StatefulPartitionedCall"^dense_126/StatefulPartitionedCall"^dense_127/StatefulPartitionedCall"^dense_128/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2`
.batch_normalization_84/StatefulPartitionedCall.batch_normalization_84/StatefulPartitionedCall2`
.batch_normalization_85/StatefulPartitionedCall.batch_normalization_85/StatefulPartitionedCall2F
!dense_126/StatefulPartitionedCall!dense_126/StatefulPartitionedCall2F
!dense_127/StatefulPartitionedCall!dense_127/StatefulPartitionedCall2F
!dense_128/StatefulPartitionedCall!dense_128/StatefulPartitionedCall:&"
 
_user_specified_name161550:&"
 
_user_specified_name161548:&"
 
_user_specified_name161534:&"
 
_user_specified_name161532:&
"
 
_user_specified_name161530:&	"
 
_user_specified_name161528:&"
 
_user_specified_name161525:&"
 
_user_specified_name161523:&"
 
_user_specified_name161510:&"
 
_user_specified_name161508:&"
 
_user_specified_name161506:&"
 
_user_specified_name161504:&"
 
_user_specified_name161501:&"
 
_user_specified_name161499:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
*__inference_dense_127_layer_call_fn_161828

inputs
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_127_layer_call_and_return_conditional_losses_161522o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161824:&"
 
_user_specified_name161822:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�%
�
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161440

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 Z
AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
AssignMovingAvg/ReadVariableOpReadVariableOp'assignmovingavg_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvgAssignSubVariableOp'assignmovingavg_readvariableop_resourceAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0\
AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
 AssignMovingAvg_1/ReadVariableOpReadVariableOp)assignmovingavg_1_readvariableop_resource*
_output_shapes
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:q
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12"
AssignMovingAvgAssignMovingAvg2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161819

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
dtype0T
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:t
batchnorm/addAddV2Cast_1/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
7__inference_batch_normalization_85_layer_call_fn_161864

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161460o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161860:&"
 
_user_specified_name161858:&"
 
_user_specified_name161856:&"
 
_user_specified_name161854:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_126_layer_call_and_return_conditional_losses_161739

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�	
�
7__inference_batch_normalization_84_layer_call_fn_161752

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *[
fVRT
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161360o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161748:&"
 
_user_specified_name161746:&"
 
_user_specified_name161744:&"
 
_user_specified_name161742:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_127_layer_call_and_return_conditional_losses_161838

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
/__inference_dense_model_42_layer_call_fn_161657
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161591o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161653:&"
 
_user_specified_name161651:&"
 
_user_specified_name161649:&"
 
_user_specified_name161647:&
"
 
_user_specified_name161645:&	"
 
_user_specified_name161643:&"
 
_user_specified_name161641:&"
 
_user_specified_name161639:&"
 
_user_specified_name161637:&"
 
_user_specified_name161635:&"
 
_user_specified_name161633:&"
 
_user_specified_name161631:&"
 
_user_specified_name161629:&"
 
_user_specified_name161627:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
*__inference_dense_128_layer_call_fn_161927

inputs
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_128_layer_call_and_return_conditional_losses_161547o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name161923:&"
 
_user_specified_name161921:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
<
input_11
serving_default_input_1:0����������<
output_10
StatefulPartitionedCall:0���������tensorflow/serving/predict:ڢ
�
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature

dense1
		norm1


dense2
	norm2
out
	optimizer

signatures"
_tf_keras_model
�
0
1
2
3
4
5
6
7
8
9
10
11
12
13"
trackable_list_wrapper
f
0
1
2
3
4
5
6
7
8
9"
trackable_list_wrapper
 "
trackable_list_wrapper
�
non_trainable_variables

layers
metrics
 layer_regularization_losses
!layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�
"trace_0
#trace_12�
/__inference_dense_model_42_layer_call_fn_161624
/__inference_dense_model_42_layer_call_fn_161657�
���
FullArgSpec
args�
jx
varargs
 
varkw
 
defaults
 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 z"trace_0z#trace_1
�
$trace_0
%trace_12�
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161554
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161591�
���
FullArgSpec
args�
jx
varargs
 
varkw
 
defaults
 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 z$trace_0z%trace_1
�B�
!__inference__wrapped_model_161326input_1"�
���
FullArgSpec
args� 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�
&	variables
'trainable_variables
(regularization_losses
)	keras_api
*__call__
*+&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
,	variables
-trainable_variables
.regularization_losses
/	keras_api
0__call__
*1&call_and_return_all_conditional_losses
2axis
	gamma
beta
moving_mean
moving_variance"
_tf_keras_layer
�
3	variables
4trainable_variables
5regularization_losses
6	keras_api
7__call__
*8&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
9	variables
:trainable_variables
;regularization_losses
<	keras_api
=__call__
*>&call_and_return_all_conditional_losses
?axis
	gamma
beta
moving_mean
moving_variance"
_tf_keras_layer
�
@	variables
Atrainable_variables
Bregularization_losses
C	keras_api
D__call__
*E&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
F
_variables
G_iterations
H_learning_rate
I_index_dict
J
_momentums
K_velocities
L_update_step_xla"
experimentalOptimizer
,
Mserving_default"
signature_map
2:0	�2dense_model_42/dense_126/kernel
+:)2dense_model_42/dense_126/bias
9:72+dense_model_42/batch_normalization_84/gamma
8:62*dense_model_42/batch_normalization_84/beta
A:? (21dense_model_42/batch_normalization_84/moving_mean
E:C (25dense_model_42/batch_normalization_84/moving_variance
1:/2dense_model_42/dense_127/kernel
+:)2dense_model_42/dense_127/bias
9:72+dense_model_42/batch_normalization_85/gamma
8:62*dense_model_42/batch_normalization_85/beta
A:? (21dense_model_42/batch_normalization_85/moving_mean
E:C (25dense_model_42/batch_normalization_85/moving_variance
1:/2dense_model_42/dense_128/kernel
+:)2dense_model_42/dense_128/bias
<
0
1
2
3"
trackable_list_wrapper
C
0
	1

2
3
4"
trackable_list_wrapper
.
N0
O1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
/__inference_dense_model_42_layer_call_fn_161624input_1"�
���
FullArgSpec
args�
jx
varargs
 
varkw
 
defaults
 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�B�
/__inference_dense_model_42_layer_call_fn_161657input_1"�
���
FullArgSpec
args�
jx
varargs
 
varkw
 
defaults
 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�B�
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161554input_1"�
���
FullArgSpec
args�
jx
varargs
 
varkw
 
defaults
 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�B�
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161591input_1"�
���
FullArgSpec
args�
jx
varargs
 
varkw
 
defaults
 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
Pnon_trainable_variables

Qlayers
Rmetrics
Slayer_regularization_losses
Tlayer_metrics
&	variables
'trainable_variables
(regularization_losses
*__call__
*+&call_and_return_all_conditional_losses
&+"call_and_return_conditional_losses"
_generic_user_object
�
Utrace_02�
*__inference_dense_126_layer_call_fn_161729�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zUtrace_0
�
Vtrace_02�
E__inference_dense_126_layer_call_and_return_conditional_losses_161739�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zVtrace_0
<
0
1
2
3"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
Wnon_trainable_variables

Xlayers
Ymetrics
Zlayer_regularization_losses
[layer_metrics
,	variables
-trainable_variables
.regularization_losses
0__call__
*1&call_and_return_all_conditional_losses
&1"call_and_return_conditional_losses"
_generic_user_object
�
\trace_0
]trace_12�
7__inference_batch_normalization_84_layer_call_fn_161752
7__inference_batch_normalization_84_layer_call_fn_161765�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z\trace_0z]trace_1
�
^trace_0
_trace_12�
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161799
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161819�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z^trace_0z_trace_1
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
`non_trainable_variables

alayers
bmetrics
clayer_regularization_losses
dlayer_metrics
3	variables
4trainable_variables
5regularization_losses
7__call__
*8&call_and_return_all_conditional_losses
&8"call_and_return_conditional_losses"
_generic_user_object
�
etrace_02�
*__inference_dense_127_layer_call_fn_161828�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zetrace_0
�
ftrace_02�
E__inference_dense_127_layer_call_and_return_conditional_losses_161838�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zftrace_0
<
0
1
2
3"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
gnon_trainable_variables

hlayers
imetrics
jlayer_regularization_losses
klayer_metrics
9	variables
:trainable_variables
;regularization_losses
=__call__
*>&call_and_return_all_conditional_losses
&>"call_and_return_conditional_losses"
_generic_user_object
�
ltrace_0
mtrace_12�
7__inference_batch_normalization_85_layer_call_fn_161851
7__inference_batch_normalization_85_layer_call_fn_161864�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zltrace_0zmtrace_1
�
ntrace_0
otrace_12�
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161898
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161918�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zntrace_0zotrace_1
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
pnon_trainable_variables

qlayers
rmetrics
slayer_regularization_losses
tlayer_metrics
@	variables
Atrainable_variables
Bregularization_losses
D__call__
*E&call_and_return_all_conditional_losses
&E"call_and_return_conditional_losses"
_generic_user_object
�
utrace_02�
*__inference_dense_128_layer_call_fn_161927�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zutrace_0
�
vtrace_02�
E__inference_dense_128_layer_call_and_return_conditional_losses_161938�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zvtrace_0
�
G0
w1
x2
y3
z4
{5
|6
}7
~8
9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20"
trackable_list_wrapper
:	 2	iteration
: 2learning_rate
 "
trackable_dict_wrapper
k
w0
y1
{2
}3
4
�5
�6
�7
�8
�9"
trackable_list_wrapper
l
x0
z1
|2
~3
�4
�5
�6
�7
�8
�9"
trackable_list_wrapper
�2��
���
FullArgSpec*
args"�

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 0
�B�
$__inference_signature_wrapper_161720input_1"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
R
�	variables
�	keras_api

�total

�count"
_tf_keras_metric
c
�	variables
�	keras_api

�total

�count
�
_fn_kwargs"
_tf_keras_metric
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
*__inference_dense_126_layer_call_fn_161729inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
E__inference_dense_126_layer_call_and_return_conditional_losses_161739inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
7__inference_batch_normalization_84_layer_call_fn_161752inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
7__inference_batch_normalization_84_layer_call_fn_161765inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161799inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161819inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
*__inference_dense_127_layer_call_fn_161828inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
E__inference_dense_127_layer_call_and_return_conditional_losses_161838inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
7__inference_batch_normalization_85_layer_call_fn_161851inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
7__inference_batch_normalization_85_layer_call_fn_161864inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161898inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161918inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
*__inference_dense_128_layer_call_fn_161927inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
E__inference_dense_128_layer_call_and_return_conditional_losses_161938inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
7:5	�2&Adam/m/dense_model_42/dense_126/kernel
7:5	�2&Adam/v/dense_model_42/dense_126/kernel
0:.2$Adam/m/dense_model_42/dense_126/bias
0:.2$Adam/v/dense_model_42/dense_126/bias
>:<22Adam/m/dense_model_42/batch_normalization_84/gamma
>:<22Adam/v/dense_model_42/batch_normalization_84/gamma
=:;21Adam/m/dense_model_42/batch_normalization_84/beta
=:;21Adam/v/dense_model_42/batch_normalization_84/beta
6:42&Adam/m/dense_model_42/dense_127/kernel
6:42&Adam/v/dense_model_42/dense_127/kernel
0:.2$Adam/m/dense_model_42/dense_127/bias
0:.2$Adam/v/dense_model_42/dense_127/bias
>:<22Adam/m/dense_model_42/batch_normalization_85/gamma
>:<22Adam/v/dense_model_42/batch_normalization_85/gamma
=:;21Adam/m/dense_model_42/batch_normalization_85/beta
=:;21Adam/v/dense_model_42/batch_normalization_85/beta
6:42&Adam/m/dense_model_42/dense_128/kernel
6:42&Adam/v/dense_model_42/dense_128/kernel
0:.2$Adam/m/dense_model_42/dense_128/bias
0:.2$Adam/v/dense_model_42/dense_128/bias
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper�
!__inference__wrapped_model_161326x1�.
'�$
"�
input_1����������
� "3�0
.
output_1"�
output_1����������
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161799m7�4
-�*
 �
inputs���������
p

 
� ",�)
"�
tensor_0���������
� �
R__inference_batch_normalization_84_layer_call_and_return_conditional_losses_161819m7�4
-�*
 �
inputs���������
p 

 
� ",�)
"�
tensor_0���������
� �
7__inference_batch_normalization_84_layer_call_fn_161752b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
7__inference_batch_normalization_84_layer_call_fn_161765b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161898m7�4
-�*
 �
inputs���������
p

 
� ",�)
"�
tensor_0���������
� �
R__inference_batch_normalization_85_layer_call_and_return_conditional_losses_161918m7�4
-�*
 �
inputs���������
p 

 
� ",�)
"�
tensor_0���������
� �
7__inference_batch_normalization_85_layer_call_fn_161851b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
7__inference_batch_normalization_85_layer_call_fn_161864b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
E__inference_dense_126_layer_call_and_return_conditional_losses_161739d0�-
&�#
!�
inputs����������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_126_layer_call_fn_161729Y0�-
&�#
!�
inputs����������
� "!�
unknown����������
E__inference_dense_127_layer_call_and_return_conditional_losses_161838c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_127_layer_call_fn_161828X/�,
%�"
 �
inputs���������
� "!�
unknown����������
E__inference_dense_128_layer_call_and_return_conditional_losses_161938c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_128_layer_call_fn_161927X/�,
%�"
 �
inputs���������
� "!�
unknown����������
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161554�A�>
'�$
"�
input_1����������
�

trainingp",�)
"�
tensor_0���������
� �
J__inference_dense_model_42_layer_call_and_return_conditional_losses_161591�A�>
'�$
"�
input_1����������
�

trainingp ",�)
"�
tensor_0���������
� �
/__inference_dense_model_42_layer_call_fn_161624vA�>
'�$
"�
input_1����������
�

trainingp"!�
unknown����������
/__inference_dense_model_42_layer_call_fn_161657vA�>
'�$
"�
input_1����������
�

trainingp "!�
unknown����������
$__inference_signature_wrapper_161720�<�9
� 
2�/
-
input_1"�
input_1����������"3�0
.
output_1"�
output_1���������