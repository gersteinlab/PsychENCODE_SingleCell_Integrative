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
%Adam/v/dense_model_103/dense_311/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/v/dense_model_103/dense_311/bias
�
9Adam/v/dense_model_103/dense_311/bias/Read/ReadVariableOpReadVariableOp%Adam/v/dense_model_103/dense_311/bias*
_output_shapes
:*
dtype0
�
%Adam/m/dense_model_103/dense_311/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/m/dense_model_103/dense_311/bias
�
9Adam/m/dense_model_103/dense_311/bias/Read/ReadVariableOpReadVariableOp%Adam/m/dense_model_103/dense_311/bias*
_output_shapes
:*
dtype0
�
'Adam/v/dense_model_103/dense_311/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/v/dense_model_103/dense_311/kernel
�
;Adam/v/dense_model_103/dense_311/kernel/Read/ReadVariableOpReadVariableOp'Adam/v/dense_model_103/dense_311/kernel*
_output_shapes

:*
dtype0
�
'Adam/m/dense_model_103/dense_311/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/m/dense_model_103/dense_311/kernel
�
;Adam/m/dense_model_103/dense_311/kernel/Read/ReadVariableOpReadVariableOp'Adam/m/dense_model_103/dense_311/kernel*
_output_shapes

:*
dtype0
�
3Adam/v/dense_model_103/batch_normalization_207/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/dense_model_103/batch_normalization_207/beta
�
GAdam/v/dense_model_103/batch_normalization_207/beta/Read/ReadVariableOpReadVariableOp3Adam/v/dense_model_103/batch_normalization_207/beta*
_output_shapes
:*
dtype0
�
3Adam/m/dense_model_103/batch_normalization_207/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/dense_model_103/batch_normalization_207/beta
�
GAdam/m/dense_model_103/batch_normalization_207/beta/Read/ReadVariableOpReadVariableOp3Adam/m/dense_model_103/batch_normalization_207/beta*
_output_shapes
:*
dtype0
�
4Adam/v/dense_model_103/batch_normalization_207/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/v/dense_model_103/batch_normalization_207/gamma
�
HAdam/v/dense_model_103/batch_normalization_207/gamma/Read/ReadVariableOpReadVariableOp4Adam/v/dense_model_103/batch_normalization_207/gamma*
_output_shapes
:*
dtype0
�
4Adam/m/dense_model_103/batch_normalization_207/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/m/dense_model_103/batch_normalization_207/gamma
�
HAdam/m/dense_model_103/batch_normalization_207/gamma/Read/ReadVariableOpReadVariableOp4Adam/m/dense_model_103/batch_normalization_207/gamma*
_output_shapes
:*
dtype0
�
%Adam/v/dense_model_103/dense_310/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/v/dense_model_103/dense_310/bias
�
9Adam/v/dense_model_103/dense_310/bias/Read/ReadVariableOpReadVariableOp%Adam/v/dense_model_103/dense_310/bias*
_output_shapes
:*
dtype0
�
%Adam/m/dense_model_103/dense_310/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/m/dense_model_103/dense_310/bias
�
9Adam/m/dense_model_103/dense_310/bias/Read/ReadVariableOpReadVariableOp%Adam/m/dense_model_103/dense_310/bias*
_output_shapes
:*
dtype0
�
'Adam/v/dense_model_103/dense_310/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/v/dense_model_103/dense_310/kernel
�
;Adam/v/dense_model_103/dense_310/kernel/Read/ReadVariableOpReadVariableOp'Adam/v/dense_model_103/dense_310/kernel*
_output_shapes

:*
dtype0
�
'Adam/m/dense_model_103/dense_310/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/m/dense_model_103/dense_310/kernel
�
;Adam/m/dense_model_103/dense_310/kernel/Read/ReadVariableOpReadVariableOp'Adam/m/dense_model_103/dense_310/kernel*
_output_shapes

:*
dtype0
�
3Adam/v/dense_model_103/batch_normalization_206/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/dense_model_103/batch_normalization_206/beta
�
GAdam/v/dense_model_103/batch_normalization_206/beta/Read/ReadVariableOpReadVariableOp3Adam/v/dense_model_103/batch_normalization_206/beta*
_output_shapes
:*
dtype0
�
3Adam/m/dense_model_103/batch_normalization_206/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/dense_model_103/batch_normalization_206/beta
�
GAdam/m/dense_model_103/batch_normalization_206/beta/Read/ReadVariableOpReadVariableOp3Adam/m/dense_model_103/batch_normalization_206/beta*
_output_shapes
:*
dtype0
�
4Adam/v/dense_model_103/batch_normalization_206/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/v/dense_model_103/batch_normalization_206/gamma
�
HAdam/v/dense_model_103/batch_normalization_206/gamma/Read/ReadVariableOpReadVariableOp4Adam/v/dense_model_103/batch_normalization_206/gamma*
_output_shapes
:*
dtype0
�
4Adam/m/dense_model_103/batch_normalization_206/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/m/dense_model_103/batch_normalization_206/gamma
�
HAdam/m/dense_model_103/batch_normalization_206/gamma/Read/ReadVariableOpReadVariableOp4Adam/m/dense_model_103/batch_normalization_206/gamma*
_output_shapes
:*
dtype0
�
%Adam/v/dense_model_103/dense_309/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/v/dense_model_103/dense_309/bias
�
9Adam/v/dense_model_103/dense_309/bias/Read/ReadVariableOpReadVariableOp%Adam/v/dense_model_103/dense_309/bias*
_output_shapes
:*
dtype0
�
%Adam/m/dense_model_103/dense_309/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/m/dense_model_103/dense_309/bias
�
9Adam/m/dense_model_103/dense_309/bias/Read/ReadVariableOpReadVariableOp%Adam/m/dense_model_103/dense_309/bias*
_output_shapes
:*
dtype0
�
'Adam/v/dense_model_103/dense_309/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*8
shared_name)'Adam/v/dense_model_103/dense_309/kernel
�
;Adam/v/dense_model_103/dense_309/kernel/Read/ReadVariableOpReadVariableOp'Adam/v/dense_model_103/dense_309/kernel*
_output_shapes
:	�*
dtype0
�
'Adam/m/dense_model_103/dense_309/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*8
shared_name)'Adam/m/dense_model_103/dense_309/kernel
�
;Adam/m/dense_model_103/dense_309/kernel/Read/ReadVariableOpReadVariableOp'Adam/m/dense_model_103/dense_309/kernel*
_output_shapes
:	�*
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
dense_model_103/dense_311/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name dense_model_103/dense_311/bias
�
2dense_model_103/dense_311/bias/Read/ReadVariableOpReadVariableOpdense_model_103/dense_311/bias*
_output_shapes
:*
dtype0
�
 dense_model_103/dense_311/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*1
shared_name" dense_model_103/dense_311/kernel
�
4dense_model_103/dense_311/kernel/Read/ReadVariableOpReadVariableOp dense_model_103/dense_311/kernel*
_output_shapes

:*
dtype0
�
7dense_model_103/batch_normalization_207/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*H
shared_name97dense_model_103/batch_normalization_207/moving_variance
�
Kdense_model_103/batch_normalization_207/moving_variance/Read/ReadVariableOpReadVariableOp7dense_model_103/batch_normalization_207/moving_variance*
_output_shapes
:*
dtype0
�
3dense_model_103/batch_normalization_207/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53dense_model_103/batch_normalization_207/moving_mean
�
Gdense_model_103/batch_normalization_207/moving_mean/Read/ReadVariableOpReadVariableOp3dense_model_103/batch_normalization_207/moving_mean*
_output_shapes
:*
dtype0
�
,dense_model_103/batch_normalization_207/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,dense_model_103/batch_normalization_207/beta
�
@dense_model_103/batch_normalization_207/beta/Read/ReadVariableOpReadVariableOp,dense_model_103/batch_normalization_207/beta*
_output_shapes
:*
dtype0
�
-dense_model_103/batch_normalization_207/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*>
shared_name/-dense_model_103/batch_normalization_207/gamma
�
Adense_model_103/batch_normalization_207/gamma/Read/ReadVariableOpReadVariableOp-dense_model_103/batch_normalization_207/gamma*
_output_shapes
:*
dtype0
�
dense_model_103/dense_310/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name dense_model_103/dense_310/bias
�
2dense_model_103/dense_310/bias/Read/ReadVariableOpReadVariableOpdense_model_103/dense_310/bias*
_output_shapes
:*
dtype0
�
 dense_model_103/dense_310/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*1
shared_name" dense_model_103/dense_310/kernel
�
4dense_model_103/dense_310/kernel/Read/ReadVariableOpReadVariableOp dense_model_103/dense_310/kernel*
_output_shapes

:*
dtype0
�
7dense_model_103/batch_normalization_206/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*H
shared_name97dense_model_103/batch_normalization_206/moving_variance
�
Kdense_model_103/batch_normalization_206/moving_variance/Read/ReadVariableOpReadVariableOp7dense_model_103/batch_normalization_206/moving_variance*
_output_shapes
:*
dtype0
�
3dense_model_103/batch_normalization_206/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53dense_model_103/batch_normalization_206/moving_mean
�
Gdense_model_103/batch_normalization_206/moving_mean/Read/ReadVariableOpReadVariableOp3dense_model_103/batch_normalization_206/moving_mean*
_output_shapes
:*
dtype0
�
,dense_model_103/batch_normalization_206/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,dense_model_103/batch_normalization_206/beta
�
@dense_model_103/batch_normalization_206/beta/Read/ReadVariableOpReadVariableOp,dense_model_103/batch_normalization_206/beta*
_output_shapes
:*
dtype0
�
-dense_model_103/batch_normalization_206/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*>
shared_name/-dense_model_103/batch_normalization_206/gamma
�
Adense_model_103/batch_normalization_206/gamma/Read/ReadVariableOpReadVariableOp-dense_model_103/batch_normalization_206/gamma*
_output_shapes
:*
dtype0
�
dense_model_103/dense_309/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name dense_model_103/dense_309/bias
�
2dense_model_103/dense_309/bias/Read/ReadVariableOpReadVariableOpdense_model_103/dense_309/bias*
_output_shapes
:*
dtype0
�
 dense_model_103/dense_309/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*1
shared_name" dense_model_103/dense_309/kernel
�
4dense_model_103/dense_309/kernel/Read/ReadVariableOpReadVariableOp dense_model_103/dense_309/kernel*
_output_shapes
:	�*
dtype0
|
serving_default_input_1Placeholder*(
_output_shapes
:����������*
dtype0*
shape:����������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1 dense_model_103/dense_309/kerneldense_model_103/dense_309/bias3dense_model_103/batch_normalization_206/moving_mean7dense_model_103/batch_normalization_206/moving_variance,dense_model_103/batch_normalization_206/beta-dense_model_103/batch_normalization_206/gamma dense_model_103/dense_310/kerneldense_model_103/dense_310/bias3dense_model_103/batch_normalization_207/moving_mean7dense_model_103/batch_normalization_207/moving_variance,dense_model_103/batch_normalization_207/beta-dense_model_103/batch_normalization_207/gamma dense_model_103/dense_311/kerneldense_model_103/dense_311/bias*
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
$__inference_signature_wrapper_390655

NoOpNoOp
�G
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�G
value�GB�G B�G
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
`Z
VARIABLE_VALUE dense_model_103/dense_309/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEdense_model_103/dense_309/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
mg
VARIABLE_VALUE-dense_model_103/batch_normalization_206/gamma&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE,dense_model_103/batch_normalization_206/beta&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE3dense_model_103/batch_normalization_206/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
wq
VARIABLE_VALUE7dense_model_103/batch_normalization_206/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUE dense_model_103/dense_310/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEdense_model_103/dense_310/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE*
mg
VARIABLE_VALUE-dense_model_103/batch_normalization_207/gamma&variables/8/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE,dense_model_103/batch_normalization_207/beta&variables/9/.ATTRIBUTES/VARIABLE_VALUE*
tn
VARIABLE_VALUE3dense_model_103/batch_normalization_207/moving_mean'variables/10/.ATTRIBUTES/VARIABLE_VALUE*
xr
VARIABLE_VALUE7dense_model_103/batch_normalization_207/moving_variance'variables/11/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUE dense_model_103/dense_311/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEdense_model_103/dense_311/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE*
 
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
rl
VARIABLE_VALUE'Adam/m/dense_model_103/dense_309/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE'Adam/v/dense_model_103/dense_309/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE%Adam/m/dense_model_103/dense_309/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE%Adam/v/dense_model_103/dense_309/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE4Adam/m/dense_model_103/batch_normalization_206/gamma1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE4Adam/v/dense_model_103/batch_normalization_206/gamma1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE3Adam/m/dense_model_103/batch_normalization_206/beta1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE3Adam/v/dense_model_103/batch_normalization_206/beta1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE'Adam/m/dense_model_103/dense_310/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/v/dense_model_103/dense_310/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/m/dense_model_103/dense_310/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/v/dense_model_103/dense_310/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUE4Adam/m/dense_model_103/batch_normalization_207/gamma2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUE4Adam/v/dense_model_103/batch_normalization_207/gamma2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/m/dense_model_103/batch_normalization_207/beta2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/v/dense_model_103/batch_normalization_207/beta2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/m/dense_model_103/dense_311/kernel2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/v/dense_model_103/dense_311/kernel2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/m/dense_model_103/dense_311/bias2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/v/dense_model_103/dense_311/bias2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUE*
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
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename dense_model_103/dense_309/kerneldense_model_103/dense_309/bias-dense_model_103/batch_normalization_206/gamma,dense_model_103/batch_normalization_206/beta3dense_model_103/batch_normalization_206/moving_mean7dense_model_103/batch_normalization_206/moving_variance dense_model_103/dense_310/kerneldense_model_103/dense_310/bias-dense_model_103/batch_normalization_207/gamma,dense_model_103/batch_normalization_207/beta3dense_model_103/batch_normalization_207/moving_mean7dense_model_103/batch_normalization_207/moving_variance dense_model_103/dense_311/kerneldense_model_103/dense_311/bias	iterationlearning_rate'Adam/m/dense_model_103/dense_309/kernel'Adam/v/dense_model_103/dense_309/kernel%Adam/m/dense_model_103/dense_309/bias%Adam/v/dense_model_103/dense_309/bias4Adam/m/dense_model_103/batch_normalization_206/gamma4Adam/v/dense_model_103/batch_normalization_206/gamma3Adam/m/dense_model_103/batch_normalization_206/beta3Adam/v/dense_model_103/batch_normalization_206/beta'Adam/m/dense_model_103/dense_310/kernel'Adam/v/dense_model_103/dense_310/kernel%Adam/m/dense_model_103/dense_310/bias%Adam/v/dense_model_103/dense_310/bias4Adam/m/dense_model_103/batch_normalization_207/gamma4Adam/v/dense_model_103/batch_normalization_207/gamma3Adam/m/dense_model_103/batch_normalization_207/beta3Adam/v/dense_model_103/batch_normalization_207/beta'Adam/m/dense_model_103/dense_311/kernel'Adam/v/dense_model_103/dense_311/kernel%Adam/m/dense_model_103/dense_311/bias%Adam/v/dense_model_103/dense_311/biastotal_1count_1totalcountConst*5
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
__inference__traced_save_391135
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename dense_model_103/dense_309/kerneldense_model_103/dense_309/bias-dense_model_103/batch_normalization_206/gamma,dense_model_103/batch_normalization_206/beta3dense_model_103/batch_normalization_206/moving_mean7dense_model_103/batch_normalization_206/moving_variance dense_model_103/dense_310/kerneldense_model_103/dense_310/bias-dense_model_103/batch_normalization_207/gamma,dense_model_103/batch_normalization_207/beta3dense_model_103/batch_normalization_207/moving_mean7dense_model_103/batch_normalization_207/moving_variance dense_model_103/dense_311/kerneldense_model_103/dense_311/bias	iterationlearning_rate'Adam/m/dense_model_103/dense_309/kernel'Adam/v/dense_model_103/dense_309/kernel%Adam/m/dense_model_103/dense_309/bias%Adam/v/dense_model_103/dense_309/bias4Adam/m/dense_model_103/batch_normalization_206/gamma4Adam/v/dense_model_103/batch_normalization_206/gamma3Adam/m/dense_model_103/batch_normalization_206/beta3Adam/v/dense_model_103/batch_normalization_206/beta'Adam/m/dense_model_103/dense_310/kernel'Adam/v/dense_model_103/dense_310/kernel%Adam/m/dense_model_103/dense_310/bias%Adam/v/dense_model_103/dense_310/bias4Adam/m/dense_model_103/batch_normalization_207/gamma4Adam/v/dense_model_103/batch_normalization_207/gamma3Adam/m/dense_model_103/batch_normalization_207/beta3Adam/v/dense_model_103/batch_normalization_207/beta'Adam/m/dense_model_103/dense_311/kernel'Adam/v/dense_model_103/dense_311/kernel%Adam/m/dense_model_103/dense_311/bias%Adam/v/dense_model_103/dense_311/biastotal_1count_1totalcount*4
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
"__inference__traced_restore_391264��
�b
�
!__inference__wrapped_model_390261
input_1K
8dense_model_103_dense_309_matmul_readvariableop_resource:	�G
9dense_model_103_dense_309_biasadd_readvariableop_resource:R
Ddense_model_103_batch_normalization_206_cast_readvariableop_resource:T
Fdense_model_103_batch_normalization_206_cast_1_readvariableop_resource:T
Fdense_model_103_batch_normalization_206_cast_2_readvariableop_resource:T
Fdense_model_103_batch_normalization_206_cast_3_readvariableop_resource:J
8dense_model_103_dense_310_matmul_readvariableop_resource:G
9dense_model_103_dense_310_biasadd_readvariableop_resource:R
Ddense_model_103_batch_normalization_207_cast_readvariableop_resource:T
Fdense_model_103_batch_normalization_207_cast_1_readvariableop_resource:T
Fdense_model_103_batch_normalization_207_cast_2_readvariableop_resource:T
Fdense_model_103_batch_normalization_207_cast_3_readvariableop_resource:J
8dense_model_103_dense_311_matmul_readvariableop_resource:G
9dense_model_103_dense_311_biasadd_readvariableop_resource:
identity��;dense_model_103/batch_normalization_206/Cast/ReadVariableOp�=dense_model_103/batch_normalization_206/Cast_1/ReadVariableOp�=dense_model_103/batch_normalization_206/Cast_2/ReadVariableOp�=dense_model_103/batch_normalization_206/Cast_3/ReadVariableOp�;dense_model_103/batch_normalization_207/Cast/ReadVariableOp�=dense_model_103/batch_normalization_207/Cast_1/ReadVariableOp�=dense_model_103/batch_normalization_207/Cast_2/ReadVariableOp�=dense_model_103/batch_normalization_207/Cast_3/ReadVariableOp�0dense_model_103/dense_309/BiasAdd/ReadVariableOp�/dense_model_103/dense_309/MatMul/ReadVariableOp�0dense_model_103/dense_310/BiasAdd/ReadVariableOp�/dense_model_103/dense_310/MatMul/ReadVariableOp�0dense_model_103/dense_311/BiasAdd/ReadVariableOp�/dense_model_103/dense_311/MatMul/ReadVariableOp�
/dense_model_103/dense_309/MatMul/ReadVariableOpReadVariableOp8dense_model_103_dense_309_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
 dense_model_103/dense_309/MatMulMatMulinput_17dense_model_103/dense_309/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
0dense_model_103/dense_309/BiasAdd/ReadVariableOpReadVariableOp9dense_model_103_dense_309_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
!dense_model_103/dense_309/BiasAddBiasAdd*dense_model_103/dense_309/MatMul:product:08dense_model_103/dense_309/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
;dense_model_103/batch_normalization_206/Cast/ReadVariableOpReadVariableOpDdense_model_103_batch_normalization_206_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
=dense_model_103/batch_normalization_206/Cast_1/ReadVariableOpReadVariableOpFdense_model_103_batch_normalization_206_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
=dense_model_103/batch_normalization_206/Cast_2/ReadVariableOpReadVariableOpFdense_model_103_batch_normalization_206_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
=dense_model_103/batch_normalization_206/Cast_3/ReadVariableOpReadVariableOpFdense_model_103_batch_normalization_206_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0|
7dense_model_103/batch_normalization_206/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
5dense_model_103/batch_normalization_206/batchnorm/addAddV2Edense_model_103/batch_normalization_206/Cast_1/ReadVariableOp:value:0@dense_model_103/batch_normalization_206/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
7dense_model_103/batch_normalization_206/batchnorm/RsqrtRsqrt9dense_model_103/batch_normalization_206/batchnorm/add:z:0*
T0*
_output_shapes
:�
5dense_model_103/batch_normalization_206/batchnorm/mulMul;dense_model_103/batch_normalization_206/batchnorm/Rsqrt:y:0Edense_model_103/batch_normalization_206/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
7dense_model_103/batch_normalization_206/batchnorm/mul_1Mul*dense_model_103/dense_309/BiasAdd:output:09dense_model_103/batch_normalization_206/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
7dense_model_103/batch_normalization_206/batchnorm/mul_2MulCdense_model_103/batch_normalization_206/Cast/ReadVariableOp:value:09dense_model_103/batch_normalization_206/batchnorm/mul:z:0*
T0*
_output_shapes
:�
5dense_model_103/batch_normalization_206/batchnorm/subSubEdense_model_103/batch_normalization_206/Cast_2/ReadVariableOp:value:0;dense_model_103/batch_normalization_206/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
7dense_model_103/batch_normalization_206/batchnorm/add_1AddV2;dense_model_103/batch_normalization_206/batchnorm/mul_1:z:09dense_model_103/batch_normalization_206/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
/dense_model_103/dense_310/MatMul/ReadVariableOpReadVariableOp8dense_model_103_dense_310_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
 dense_model_103/dense_310/MatMulMatMul;dense_model_103/batch_normalization_206/batchnorm/add_1:z:07dense_model_103/dense_310/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
0dense_model_103/dense_310/BiasAdd/ReadVariableOpReadVariableOp9dense_model_103_dense_310_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
!dense_model_103/dense_310/BiasAddBiasAdd*dense_model_103/dense_310/MatMul:product:08dense_model_103/dense_310/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
;dense_model_103/batch_normalization_207/Cast/ReadVariableOpReadVariableOpDdense_model_103_batch_normalization_207_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
=dense_model_103/batch_normalization_207/Cast_1/ReadVariableOpReadVariableOpFdense_model_103_batch_normalization_207_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
=dense_model_103/batch_normalization_207/Cast_2/ReadVariableOpReadVariableOpFdense_model_103_batch_normalization_207_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
=dense_model_103/batch_normalization_207/Cast_3/ReadVariableOpReadVariableOpFdense_model_103_batch_normalization_207_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0|
7dense_model_103/batch_normalization_207/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
5dense_model_103/batch_normalization_207/batchnorm/addAddV2Edense_model_103/batch_normalization_207/Cast_1/ReadVariableOp:value:0@dense_model_103/batch_normalization_207/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
7dense_model_103/batch_normalization_207/batchnorm/RsqrtRsqrt9dense_model_103/batch_normalization_207/batchnorm/add:z:0*
T0*
_output_shapes
:�
5dense_model_103/batch_normalization_207/batchnorm/mulMul;dense_model_103/batch_normalization_207/batchnorm/Rsqrt:y:0Edense_model_103/batch_normalization_207/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
7dense_model_103/batch_normalization_207/batchnorm/mul_1Mul*dense_model_103/dense_310/BiasAdd:output:09dense_model_103/batch_normalization_207/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
7dense_model_103/batch_normalization_207/batchnorm/mul_2MulCdense_model_103/batch_normalization_207/Cast/ReadVariableOp:value:09dense_model_103/batch_normalization_207/batchnorm/mul:z:0*
T0*
_output_shapes
:�
5dense_model_103/batch_normalization_207/batchnorm/subSubEdense_model_103/batch_normalization_207/Cast_2/ReadVariableOp:value:0;dense_model_103/batch_normalization_207/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
7dense_model_103/batch_normalization_207/batchnorm/add_1AddV2;dense_model_103/batch_normalization_207/batchnorm/mul_1:z:09dense_model_103/batch_normalization_207/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
/dense_model_103/dense_311/MatMul/ReadVariableOpReadVariableOp8dense_model_103_dense_311_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
 dense_model_103/dense_311/MatMulMatMul;dense_model_103/batch_normalization_207/batchnorm/add_1:z:07dense_model_103/dense_311/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
0dense_model_103/dense_311/BiasAdd/ReadVariableOpReadVariableOp9dense_model_103_dense_311_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
!dense_model_103/dense_311/BiasAddBiasAdd*dense_model_103/dense_311/MatMul:product:08dense_model_103/dense_311/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
!dense_model_103/dense_311/SigmoidSigmoid*dense_model_103/dense_311/BiasAdd:output:0*
T0*'
_output_shapes
:���������t
IdentityIdentity%dense_model_103/dense_311/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp<^dense_model_103/batch_normalization_206/Cast/ReadVariableOp>^dense_model_103/batch_normalization_206/Cast_1/ReadVariableOp>^dense_model_103/batch_normalization_206/Cast_2/ReadVariableOp>^dense_model_103/batch_normalization_206/Cast_3/ReadVariableOp<^dense_model_103/batch_normalization_207/Cast/ReadVariableOp>^dense_model_103/batch_normalization_207/Cast_1/ReadVariableOp>^dense_model_103/batch_normalization_207/Cast_2/ReadVariableOp>^dense_model_103/batch_normalization_207/Cast_3/ReadVariableOp1^dense_model_103/dense_309/BiasAdd/ReadVariableOp0^dense_model_103/dense_309/MatMul/ReadVariableOp1^dense_model_103/dense_310/BiasAdd/ReadVariableOp0^dense_model_103/dense_310/MatMul/ReadVariableOp1^dense_model_103/dense_311/BiasAdd/ReadVariableOp0^dense_model_103/dense_311/MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2z
;dense_model_103/batch_normalization_206/Cast/ReadVariableOp;dense_model_103/batch_normalization_206/Cast/ReadVariableOp2~
=dense_model_103/batch_normalization_206/Cast_1/ReadVariableOp=dense_model_103/batch_normalization_206/Cast_1/ReadVariableOp2~
=dense_model_103/batch_normalization_206/Cast_2/ReadVariableOp=dense_model_103/batch_normalization_206/Cast_2/ReadVariableOp2~
=dense_model_103/batch_normalization_206/Cast_3/ReadVariableOp=dense_model_103/batch_normalization_206/Cast_3/ReadVariableOp2z
;dense_model_103/batch_normalization_207/Cast/ReadVariableOp;dense_model_103/batch_normalization_207/Cast/ReadVariableOp2~
=dense_model_103/batch_normalization_207/Cast_1/ReadVariableOp=dense_model_103/batch_normalization_207/Cast_1/ReadVariableOp2~
=dense_model_103/batch_normalization_207/Cast_2/ReadVariableOp=dense_model_103/batch_normalization_207/Cast_2/ReadVariableOp2~
=dense_model_103/batch_normalization_207/Cast_3/ReadVariableOp=dense_model_103/batch_normalization_207/Cast_3/ReadVariableOp2d
0dense_model_103/dense_309/BiasAdd/ReadVariableOp0dense_model_103/dense_309/BiasAdd/ReadVariableOp2b
/dense_model_103/dense_309/MatMul/ReadVariableOp/dense_model_103/dense_309/MatMul/ReadVariableOp2d
0dense_model_103/dense_310/BiasAdd/ReadVariableOp0dense_model_103/dense_310/BiasAdd/ReadVariableOp2b
/dense_model_103/dense_310/MatMul/ReadVariableOp/dense_model_103/dense_310/MatMul/ReadVariableOp2d
0dense_model_103/dense_311/BiasAdd/ReadVariableOp0dense_model_103/dense_311/BiasAdd/ReadVariableOp2b
/dense_model_103/dense_311/MatMul/ReadVariableOp/dense_model_103/dense_311/MatMul/ReadVariableOp:($
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
�	
�
8__inference_batch_normalization_207_layer_call_fn_390786

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390375o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name390782:&"
 
_user_specified_name390780:&"
 
_user_specified_name390778:&"
 
_user_specified_name390776:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
*__inference_dense_311_layer_call_fn_390862

inputs
unknown:
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
E__inference_dense_311_layer_call_and_return_conditional_losses_390482o
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
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name390858:&"
 
_user_specified_name390856:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_309_layer_call_and_return_conditional_losses_390674

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
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
0__inference_dense_model_103_layer_call_fn_390559
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

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
GPU 2J 8� *T
fORM
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390489o
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
 
_user_specified_name390555:&"
 
_user_specified_name390553:&"
 
_user_specified_name390551:&"
 
_user_specified_name390549:&
"
 
_user_specified_name390547:&	"
 
_user_specified_name390545:&"
 
_user_specified_name390543:&"
 
_user_specified_name390541:&"
 
_user_specified_name390539:&"
 
_user_specified_name390537:&"
 
_user_specified_name390535:&"
 
_user_specified_name390533:&"
 
_user_specified_name390531:&"
 
_user_specified_name390529:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
*__inference_dense_310_layer_call_fn_390763

inputs
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_310_layer_call_and_return_conditional_losses_390457o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name390759:&"
 
_user_specified_name390757:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
8__inference_batch_normalization_206_layer_call_fn_390700

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390315o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name390696:&"
 
_user_specified_name390694:&"
 
_user_specified_name390692:&"
 
_user_specified_name390690:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
E__inference_dense_311_layer_call_and_return_conditional_losses_390482

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
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
:���������: : 20
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
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390853

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
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
:���������
 
_user_specified_nameinputs
�%
�
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390833

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
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
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
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
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
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
:���������
 
_user_specified_nameinputs
�
�
0__inference_dense_model_103_layer_call_fn_390592
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:
identity��StatefulPartitionedCall�
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
GPU 2J 8� *T
fORM
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390526o
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
 
_user_specified_name390588:&"
 
_user_specified_name390586:&"
 
_user_specified_name390584:&"
 
_user_specified_name390582:&
"
 
_user_specified_name390580:&	"
 
_user_specified_name390578:&"
 
_user_specified_name390576:&"
 
_user_specified_name390574:&"
 
_user_specified_name390572:&"
 
_user_specified_name390570:&"
 
_user_specified_name390568:&"
 
_user_specified_name390566:&"
 
_user_specified_name390564:&"
 
_user_specified_name390562:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�	
�
E__inference_dense_309_layer_call_and_return_conditional_losses_390433

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
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
8__inference_batch_normalization_207_layer_call_fn_390799

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390395o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name390795:&"
 
_user_specified_name390793:&"
 
_user_specified_name390791:&"
 
_user_specified_name390789:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_310_layer_call_and_return_conditional_losses_390773

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
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
:���������
 
_user_specified_nameinputs
�%
�
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390734

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
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
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
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
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
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
:���������
 
_user_specified_nameinputs
�
�
*__inference_dense_309_layer_call_fn_390664

inputs
unknown:	�
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_309_layer_call_and_return_conditional_losses_390433o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name390660:&"
 
_user_specified_name390658:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�	
�
E__inference_dense_310_layer_call_and_return_conditional_losses_390457

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
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
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390315

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
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
:���������
 
_user_specified_nameinputs
��
�*
__inference__traced_save_391135
file_prefixJ
7read_disablecopyonread_dense_model_103_dense_309_kernel:	�E
7read_1_disablecopyonread_dense_model_103_dense_309_bias:T
Fread_2_disablecopyonread_dense_model_103_batch_normalization_206_gamma:S
Eread_3_disablecopyonread_dense_model_103_batch_normalization_206_beta:Z
Lread_4_disablecopyonread_dense_model_103_batch_normalization_206_moving_mean:^
Pread_5_disablecopyonread_dense_model_103_batch_normalization_206_moving_variance:K
9read_6_disablecopyonread_dense_model_103_dense_310_kernel:E
7read_7_disablecopyonread_dense_model_103_dense_310_bias:T
Fread_8_disablecopyonread_dense_model_103_batch_normalization_207_gamma:S
Eread_9_disablecopyonread_dense_model_103_batch_normalization_207_beta:[
Mread_10_disablecopyonread_dense_model_103_batch_normalization_207_moving_mean:_
Qread_11_disablecopyonread_dense_model_103_batch_normalization_207_moving_variance:L
:read_12_disablecopyonread_dense_model_103_dense_311_kernel:F
8read_13_disablecopyonread_dense_model_103_dense_311_bias:-
#read_14_disablecopyonread_iteration:	 1
'read_15_disablecopyonread_learning_rate: T
Aread_16_disablecopyonread_adam_m_dense_model_103_dense_309_kernel:	�T
Aread_17_disablecopyonread_adam_v_dense_model_103_dense_309_kernel:	�M
?read_18_disablecopyonread_adam_m_dense_model_103_dense_309_bias:M
?read_19_disablecopyonread_adam_v_dense_model_103_dense_309_bias:\
Nread_20_disablecopyonread_adam_m_dense_model_103_batch_normalization_206_gamma:\
Nread_21_disablecopyonread_adam_v_dense_model_103_batch_normalization_206_gamma:[
Mread_22_disablecopyonread_adam_m_dense_model_103_batch_normalization_206_beta:[
Mread_23_disablecopyonread_adam_v_dense_model_103_batch_normalization_206_beta:S
Aread_24_disablecopyonread_adam_m_dense_model_103_dense_310_kernel:S
Aread_25_disablecopyonread_adam_v_dense_model_103_dense_310_kernel:M
?read_26_disablecopyonread_adam_m_dense_model_103_dense_310_bias:M
?read_27_disablecopyonread_adam_v_dense_model_103_dense_310_bias:\
Nread_28_disablecopyonread_adam_m_dense_model_103_batch_normalization_207_gamma:\
Nread_29_disablecopyonread_adam_v_dense_model_103_batch_normalization_207_gamma:[
Mread_30_disablecopyonread_adam_m_dense_model_103_batch_normalization_207_beta:[
Mread_31_disablecopyonread_adam_v_dense_model_103_batch_normalization_207_beta:S
Aread_32_disablecopyonread_adam_m_dense_model_103_dense_311_kernel:S
Aread_33_disablecopyonread_adam_v_dense_model_103_dense_311_kernel:M
?read_34_disablecopyonread_adam_m_dense_model_103_dense_311_bias:M
?read_35_disablecopyonread_adam_v_dense_model_103_dense_311_bias:+
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
Read/DisableCopyOnReadDisableCopyOnRead7read_disablecopyonread_dense_model_103_dense_309_kernel"/device:CPU:0*
_output_shapes
 �
Read/ReadVariableOpReadVariableOp7read_disablecopyonread_dense_model_103_dense_309_kernel^Read/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0j
IdentityIdentityRead/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�b

Identity_1IdentityIdentity:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_1/DisableCopyOnReadDisableCopyOnRead7read_1_disablecopyonread_dense_model_103_dense_309_bias"/device:CPU:0*
_output_shapes
 �
Read_1/ReadVariableOpReadVariableOp7read_1_disablecopyonread_dense_model_103_dense_309_bias^Read_1/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_2IdentityRead_1/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_3IdentityIdentity_2:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_2/DisableCopyOnReadDisableCopyOnReadFread_2_disablecopyonread_dense_model_103_batch_normalization_206_gamma"/device:CPU:0*
_output_shapes
 �
Read_2/ReadVariableOpReadVariableOpFread_2_disablecopyonread_dense_model_103_batch_normalization_206_gamma^Read_2/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_4IdentityRead_2/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_5IdentityIdentity_4:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_3/DisableCopyOnReadDisableCopyOnReadEread_3_disablecopyonread_dense_model_103_batch_normalization_206_beta"/device:CPU:0*
_output_shapes
 �
Read_3/ReadVariableOpReadVariableOpEread_3_disablecopyonread_dense_model_103_batch_normalization_206_beta^Read_3/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_6IdentityRead_3/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_7IdentityIdentity_6:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_4/DisableCopyOnReadDisableCopyOnReadLread_4_disablecopyonread_dense_model_103_batch_normalization_206_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_4/ReadVariableOpReadVariableOpLread_4_disablecopyonread_dense_model_103_batch_normalization_206_moving_mean^Read_4/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_8IdentityRead_4/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_9IdentityIdentity_8:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_5/DisableCopyOnReadDisableCopyOnReadPread_5_disablecopyonread_dense_model_103_batch_normalization_206_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_5/ReadVariableOpReadVariableOpPread_5_disablecopyonread_dense_model_103_batch_normalization_206_moving_variance^Read_5/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_10IdentityRead_5/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_11IdentityIdentity_10:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_6/DisableCopyOnReadDisableCopyOnRead9read_6_disablecopyonread_dense_model_103_dense_310_kernel"/device:CPU:0*
_output_shapes
 �
Read_6/ReadVariableOpReadVariableOp9read_6_disablecopyonread_dense_model_103_dense_310_kernel^Read_6/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0n
Identity_12IdentityRead_6/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_13IdentityIdentity_12:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_7/DisableCopyOnReadDisableCopyOnRead7read_7_disablecopyonread_dense_model_103_dense_310_bias"/device:CPU:0*
_output_shapes
 �
Read_7/ReadVariableOpReadVariableOp7read_7_disablecopyonread_dense_model_103_dense_310_bias^Read_7/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_14IdentityRead_7/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_15IdentityIdentity_14:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_8/DisableCopyOnReadDisableCopyOnReadFread_8_disablecopyonread_dense_model_103_batch_normalization_207_gamma"/device:CPU:0*
_output_shapes
 �
Read_8/ReadVariableOpReadVariableOpFread_8_disablecopyonread_dense_model_103_batch_normalization_207_gamma^Read_8/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_16IdentityRead_8/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_17IdentityIdentity_16:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_9/DisableCopyOnReadDisableCopyOnReadEread_9_disablecopyonread_dense_model_103_batch_normalization_207_beta"/device:CPU:0*
_output_shapes
 �
Read_9/ReadVariableOpReadVariableOpEread_9_disablecopyonread_dense_model_103_batch_normalization_207_beta^Read_9/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_18IdentityRead_9/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_19IdentityIdentity_18:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_10/DisableCopyOnReadDisableCopyOnReadMread_10_disablecopyonread_dense_model_103_batch_normalization_207_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_10/ReadVariableOpReadVariableOpMread_10_disablecopyonread_dense_model_103_batch_normalization_207_moving_mean^Read_10/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_20IdentityRead_10/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_21IdentityIdentity_20:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_11/DisableCopyOnReadDisableCopyOnReadQread_11_disablecopyonread_dense_model_103_batch_normalization_207_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_11/ReadVariableOpReadVariableOpQread_11_disablecopyonread_dense_model_103_batch_normalization_207_moving_variance^Read_11/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_22IdentityRead_11/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_23IdentityIdentity_22:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_12/DisableCopyOnReadDisableCopyOnRead:read_12_disablecopyonread_dense_model_103_dense_311_kernel"/device:CPU:0*
_output_shapes
 �
Read_12/ReadVariableOpReadVariableOp:read_12_disablecopyonread_dense_model_103_dense_311_kernel^Read_12/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_24IdentityRead_12/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_25IdentityIdentity_24:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_13/DisableCopyOnReadDisableCopyOnRead8read_13_disablecopyonread_dense_model_103_dense_311_bias"/device:CPU:0*
_output_shapes
 �
Read_13/ReadVariableOpReadVariableOp8read_13_disablecopyonread_dense_model_103_dense_311_bias^Read_13/DisableCopyOnRead"/device:CPU:0*
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
Read_16/DisableCopyOnReadDisableCopyOnReadAread_16_disablecopyonread_adam_m_dense_model_103_dense_309_kernel"/device:CPU:0*
_output_shapes
 �
Read_16/ReadVariableOpReadVariableOpAread_16_disablecopyonread_adam_m_dense_model_103_dense_309_kernel^Read_16/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0p
Identity_32IdentityRead_16/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�f
Identity_33IdentityIdentity_32:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_17/DisableCopyOnReadDisableCopyOnReadAread_17_disablecopyonread_adam_v_dense_model_103_dense_309_kernel"/device:CPU:0*
_output_shapes
 �
Read_17/ReadVariableOpReadVariableOpAread_17_disablecopyonread_adam_v_dense_model_103_dense_309_kernel^Read_17/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0p
Identity_34IdentityRead_17/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�f
Identity_35IdentityIdentity_34:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_18/DisableCopyOnReadDisableCopyOnRead?read_18_disablecopyonread_adam_m_dense_model_103_dense_309_bias"/device:CPU:0*
_output_shapes
 �
Read_18/ReadVariableOpReadVariableOp?read_18_disablecopyonread_adam_m_dense_model_103_dense_309_bias^Read_18/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_36IdentityRead_18/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_37IdentityIdentity_36:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_19/DisableCopyOnReadDisableCopyOnRead?read_19_disablecopyonread_adam_v_dense_model_103_dense_309_bias"/device:CPU:0*
_output_shapes
 �
Read_19/ReadVariableOpReadVariableOp?read_19_disablecopyonread_adam_v_dense_model_103_dense_309_bias^Read_19/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_38IdentityRead_19/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_39IdentityIdentity_38:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_20/DisableCopyOnReadDisableCopyOnReadNread_20_disablecopyonread_adam_m_dense_model_103_batch_normalization_206_gamma"/device:CPU:0*
_output_shapes
 �
Read_20/ReadVariableOpReadVariableOpNread_20_disablecopyonread_adam_m_dense_model_103_batch_normalization_206_gamma^Read_20/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_40IdentityRead_20/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_41IdentityIdentity_40:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_21/DisableCopyOnReadDisableCopyOnReadNread_21_disablecopyonread_adam_v_dense_model_103_batch_normalization_206_gamma"/device:CPU:0*
_output_shapes
 �
Read_21/ReadVariableOpReadVariableOpNread_21_disablecopyonread_adam_v_dense_model_103_batch_normalization_206_gamma^Read_21/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_42IdentityRead_21/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_43IdentityIdentity_42:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_22/DisableCopyOnReadDisableCopyOnReadMread_22_disablecopyonread_adam_m_dense_model_103_batch_normalization_206_beta"/device:CPU:0*
_output_shapes
 �
Read_22/ReadVariableOpReadVariableOpMread_22_disablecopyonread_adam_m_dense_model_103_batch_normalization_206_beta^Read_22/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_44IdentityRead_22/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_45IdentityIdentity_44:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_23/DisableCopyOnReadDisableCopyOnReadMread_23_disablecopyonread_adam_v_dense_model_103_batch_normalization_206_beta"/device:CPU:0*
_output_shapes
 �
Read_23/ReadVariableOpReadVariableOpMread_23_disablecopyonread_adam_v_dense_model_103_batch_normalization_206_beta^Read_23/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_46IdentityRead_23/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_47IdentityIdentity_46:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_24/DisableCopyOnReadDisableCopyOnReadAread_24_disablecopyonread_adam_m_dense_model_103_dense_310_kernel"/device:CPU:0*
_output_shapes
 �
Read_24/ReadVariableOpReadVariableOpAread_24_disablecopyonread_adam_m_dense_model_103_dense_310_kernel^Read_24/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_48IdentityRead_24/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_49IdentityIdentity_48:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_25/DisableCopyOnReadDisableCopyOnReadAread_25_disablecopyonread_adam_v_dense_model_103_dense_310_kernel"/device:CPU:0*
_output_shapes
 �
Read_25/ReadVariableOpReadVariableOpAread_25_disablecopyonread_adam_v_dense_model_103_dense_310_kernel^Read_25/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_50IdentityRead_25/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_51IdentityIdentity_50:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_26/DisableCopyOnReadDisableCopyOnRead?read_26_disablecopyonread_adam_m_dense_model_103_dense_310_bias"/device:CPU:0*
_output_shapes
 �
Read_26/ReadVariableOpReadVariableOp?read_26_disablecopyonread_adam_m_dense_model_103_dense_310_bias^Read_26/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_52IdentityRead_26/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_53IdentityIdentity_52:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_27/DisableCopyOnReadDisableCopyOnRead?read_27_disablecopyonread_adam_v_dense_model_103_dense_310_bias"/device:CPU:0*
_output_shapes
 �
Read_27/ReadVariableOpReadVariableOp?read_27_disablecopyonread_adam_v_dense_model_103_dense_310_bias^Read_27/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_54IdentityRead_27/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_55IdentityIdentity_54:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_28/DisableCopyOnReadDisableCopyOnReadNread_28_disablecopyonread_adam_m_dense_model_103_batch_normalization_207_gamma"/device:CPU:0*
_output_shapes
 �
Read_28/ReadVariableOpReadVariableOpNread_28_disablecopyonread_adam_m_dense_model_103_batch_normalization_207_gamma^Read_28/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_56IdentityRead_28/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_57IdentityIdentity_56:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_29/DisableCopyOnReadDisableCopyOnReadNread_29_disablecopyonread_adam_v_dense_model_103_batch_normalization_207_gamma"/device:CPU:0*
_output_shapes
 �
Read_29/ReadVariableOpReadVariableOpNread_29_disablecopyonread_adam_v_dense_model_103_batch_normalization_207_gamma^Read_29/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_58IdentityRead_29/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_59IdentityIdentity_58:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_30/DisableCopyOnReadDisableCopyOnReadMread_30_disablecopyonread_adam_m_dense_model_103_batch_normalization_207_beta"/device:CPU:0*
_output_shapes
 �
Read_30/ReadVariableOpReadVariableOpMread_30_disablecopyonread_adam_m_dense_model_103_batch_normalization_207_beta^Read_30/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_60IdentityRead_30/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_61IdentityIdentity_60:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_31/DisableCopyOnReadDisableCopyOnReadMread_31_disablecopyonread_adam_v_dense_model_103_batch_normalization_207_beta"/device:CPU:0*
_output_shapes
 �
Read_31/ReadVariableOpReadVariableOpMread_31_disablecopyonread_adam_v_dense_model_103_batch_normalization_207_beta^Read_31/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_62IdentityRead_31/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_63IdentityIdentity_62:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_32/DisableCopyOnReadDisableCopyOnReadAread_32_disablecopyonread_adam_m_dense_model_103_dense_311_kernel"/device:CPU:0*
_output_shapes
 �
Read_32/ReadVariableOpReadVariableOpAread_32_disablecopyonread_adam_m_dense_model_103_dense_311_kernel^Read_32/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_64IdentityRead_32/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_65IdentityIdentity_64:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_33/DisableCopyOnReadDisableCopyOnReadAread_33_disablecopyonread_adam_v_dense_model_103_dense_311_kernel"/device:CPU:0*
_output_shapes
 �
Read_33/ReadVariableOpReadVariableOpAread_33_disablecopyonread_adam_v_dense_model_103_dense_311_kernel^Read_33/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_66IdentityRead_33/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_67IdentityIdentity_66:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_34/DisableCopyOnReadDisableCopyOnRead?read_34_disablecopyonread_adam_m_dense_model_103_dense_311_bias"/device:CPU:0*
_output_shapes
 �
Read_34/ReadVariableOpReadVariableOp?read_34_disablecopyonread_adam_m_dense_model_103_dense_311_bias^Read_34/DisableCopyOnRead"/device:CPU:0*
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
Read_35/DisableCopyOnReadDisableCopyOnRead?read_35_disablecopyonread_adam_v_dense_model_103_dense_311_bias"/device:CPU:0*
_output_shapes
 �
Read_35/ReadVariableOpReadVariableOp?read_35_disablecopyonread_adam_v_dense_model_103_dense_311_bias^Read_35/DisableCopyOnRead"/device:CPU:0*
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
_user_specified_name	total_1:E$A
?
_user_specified_name'%Adam/v/dense_model_103/dense_311/bias:E#A
?
_user_specified_name'%Adam/m/dense_model_103/dense_311/bias:G"C
A
_user_specified_name)'Adam/v/dense_model_103/dense_311/kernel:G!C
A
_user_specified_name)'Adam/m/dense_model_103/dense_311/kernel:S O
M
_user_specified_name53Adam/v/dense_model_103/batch_normalization_207/beta:SO
M
_user_specified_name53Adam/m/dense_model_103/batch_normalization_207/beta:TP
N
_user_specified_name64Adam/v/dense_model_103/batch_normalization_207/gamma:TP
N
_user_specified_name64Adam/m/dense_model_103/batch_normalization_207/gamma:EA
?
_user_specified_name'%Adam/v/dense_model_103/dense_310/bias:EA
?
_user_specified_name'%Adam/m/dense_model_103/dense_310/bias:GC
A
_user_specified_name)'Adam/v/dense_model_103/dense_310/kernel:GC
A
_user_specified_name)'Adam/m/dense_model_103/dense_310/kernel:SO
M
_user_specified_name53Adam/v/dense_model_103/batch_normalization_206/beta:SO
M
_user_specified_name53Adam/m/dense_model_103/batch_normalization_206/beta:TP
N
_user_specified_name64Adam/v/dense_model_103/batch_normalization_206/gamma:TP
N
_user_specified_name64Adam/m/dense_model_103/batch_normalization_206/gamma:EA
?
_user_specified_name'%Adam/v/dense_model_103/dense_309/bias:EA
?
_user_specified_name'%Adam/m/dense_model_103/dense_309/bias:GC
A
_user_specified_name)'Adam/v/dense_model_103/dense_309/kernel:GC
A
_user_specified_name)'Adam/m/dense_model_103/dense_309/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:>:
8
_user_specified_name dense_model_103/dense_311/bias:@<
:
_user_specified_name" dense_model_103/dense_311/kernel:WS
Q
_user_specified_name97dense_model_103/batch_normalization_207/moving_variance:SO
M
_user_specified_name53dense_model_103/batch_normalization_207/moving_mean:L
H
F
_user_specified_name.,dense_model_103/batch_normalization_207/beta:M	I
G
_user_specified_name/-dense_model_103/batch_normalization_207/gamma:>:
8
_user_specified_name dense_model_103/dense_310/bias:@<
:
_user_specified_name" dense_model_103/dense_310/kernel:WS
Q
_user_specified_name97dense_model_103/batch_normalization_206/moving_variance:SO
M
_user_specified_name53dense_model_103/batch_normalization_206/moving_mean:LH
F
_user_specified_name.,dense_model_103/batch_normalization_206/beta:MI
G
_user_specified_name/-dense_model_103/batch_normalization_206/gamma:>:
8
_user_specified_name dense_model_103/dense_309/bias:@<
:
_user_specified_name" dense_model_103/dense_309/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�%
�
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390375

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
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
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
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
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
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
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390754

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
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
:���������
 
_user_specified_nameinputs
�%
�
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390526
input_1#
dense_309_390492:	�
dense_309_390494:,
batch_normalization_206_390497:,
batch_normalization_206_390499:,
batch_normalization_206_390501:,
batch_normalization_206_390503:"
dense_310_390506:
dense_310_390508:,
batch_normalization_207_390511:,
batch_normalization_207_390513:,
batch_normalization_207_390515:,
batch_normalization_207_390517:"
dense_311_390520:
dense_311_390522:
identity��/batch_normalization_206/StatefulPartitionedCall�/batch_normalization_207/StatefulPartitionedCall�!dense_309/StatefulPartitionedCall�!dense_310/StatefulPartitionedCall�!dense_311/StatefulPartitionedCall�
!dense_309/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_309_390492dense_309_390494*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_309_layer_call_and_return_conditional_losses_390433�
/batch_normalization_206/StatefulPartitionedCallStatefulPartitionedCall*dense_309/StatefulPartitionedCall:output:0batch_normalization_206_390497batch_normalization_206_390499batch_normalization_206_390501batch_normalization_206_390503*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390315�
!dense_310/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_206/StatefulPartitionedCall:output:0dense_310_390506dense_310_390508*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_310_layer_call_and_return_conditional_losses_390457�
/batch_normalization_207/StatefulPartitionedCallStatefulPartitionedCall*dense_310/StatefulPartitionedCall:output:0batch_normalization_207_390511batch_normalization_207_390513batch_normalization_207_390515batch_normalization_207_390517*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390395�
!dense_311/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_207/StatefulPartitionedCall:output:0dense_311_390520dense_311_390522*
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
E__inference_dense_311_layer_call_and_return_conditional_losses_390482y
IdentityIdentity*dense_311/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp0^batch_normalization_206/StatefulPartitionedCall0^batch_normalization_207/StatefulPartitionedCall"^dense_309/StatefulPartitionedCall"^dense_310/StatefulPartitionedCall"^dense_311/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2b
/batch_normalization_206/StatefulPartitionedCall/batch_normalization_206/StatefulPartitionedCall2b
/batch_normalization_207/StatefulPartitionedCall/batch_normalization_207/StatefulPartitionedCall2F
!dense_309/StatefulPartitionedCall!dense_309/StatefulPartitionedCall2F
!dense_310/StatefulPartitionedCall!dense_310/StatefulPartitionedCall2F
!dense_311/StatefulPartitionedCall!dense_311/StatefulPartitionedCall:&"
 
_user_specified_name390522:&"
 
_user_specified_name390520:&"
 
_user_specified_name390517:&"
 
_user_specified_name390515:&
"
 
_user_specified_name390513:&	"
 
_user_specified_name390511:&"
 
_user_specified_name390508:&"
 
_user_specified_name390506:&"
 
_user_specified_name390503:&"
 
_user_specified_name390501:&"
 
_user_specified_name390499:&"
 
_user_specified_name390497:&"
 
_user_specified_name390494:&"
 
_user_specified_name390492:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�

�
E__inference_dense_311_layer_call_and_return_conditional_losses_390873

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
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
:���������: : 20
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
:���������
 
_user_specified_nameinputs
�%
�
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390295

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
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
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
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
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
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
:���������
 
_user_specified_nameinputs
�	
�
8__inference_batch_normalization_206_layer_call_fn_390687

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390295o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name390683:&"
 
_user_specified_name390681:&"
 
_user_specified_name390679:&"
 
_user_specified_name390677:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
��
�
"__inference__traced_restore_391264
file_prefixD
1assignvariableop_dense_model_103_dense_309_kernel:	�?
1assignvariableop_1_dense_model_103_dense_309_bias:N
@assignvariableop_2_dense_model_103_batch_normalization_206_gamma:M
?assignvariableop_3_dense_model_103_batch_normalization_206_beta:T
Fassignvariableop_4_dense_model_103_batch_normalization_206_moving_mean:X
Jassignvariableop_5_dense_model_103_batch_normalization_206_moving_variance:E
3assignvariableop_6_dense_model_103_dense_310_kernel:?
1assignvariableop_7_dense_model_103_dense_310_bias:N
@assignvariableop_8_dense_model_103_batch_normalization_207_gamma:M
?assignvariableop_9_dense_model_103_batch_normalization_207_beta:U
Gassignvariableop_10_dense_model_103_batch_normalization_207_moving_mean:Y
Kassignvariableop_11_dense_model_103_batch_normalization_207_moving_variance:F
4assignvariableop_12_dense_model_103_dense_311_kernel:@
2assignvariableop_13_dense_model_103_dense_311_bias:'
assignvariableop_14_iteration:	 +
!assignvariableop_15_learning_rate: N
;assignvariableop_16_adam_m_dense_model_103_dense_309_kernel:	�N
;assignvariableop_17_adam_v_dense_model_103_dense_309_kernel:	�G
9assignvariableop_18_adam_m_dense_model_103_dense_309_bias:G
9assignvariableop_19_adam_v_dense_model_103_dense_309_bias:V
Hassignvariableop_20_adam_m_dense_model_103_batch_normalization_206_gamma:V
Hassignvariableop_21_adam_v_dense_model_103_batch_normalization_206_gamma:U
Gassignvariableop_22_adam_m_dense_model_103_batch_normalization_206_beta:U
Gassignvariableop_23_adam_v_dense_model_103_batch_normalization_206_beta:M
;assignvariableop_24_adam_m_dense_model_103_dense_310_kernel:M
;assignvariableop_25_adam_v_dense_model_103_dense_310_kernel:G
9assignvariableop_26_adam_m_dense_model_103_dense_310_bias:G
9assignvariableop_27_adam_v_dense_model_103_dense_310_bias:V
Hassignvariableop_28_adam_m_dense_model_103_batch_normalization_207_gamma:V
Hassignvariableop_29_adam_v_dense_model_103_batch_normalization_207_gamma:U
Gassignvariableop_30_adam_m_dense_model_103_batch_normalization_207_beta:U
Gassignvariableop_31_adam_v_dense_model_103_batch_normalization_207_beta:M
;assignvariableop_32_adam_m_dense_model_103_dense_311_kernel:M
;assignvariableop_33_adam_v_dense_model_103_dense_311_kernel:G
9assignvariableop_34_adam_m_dense_model_103_dense_311_bias:G
9assignvariableop_35_adam_v_dense_model_103_dense_311_bias:%
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
AssignVariableOpAssignVariableOp1assignvariableop_dense_model_103_dense_309_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp1assignvariableop_1_dense_model_103_dense_309_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOp@assignvariableop_2_dense_model_103_batch_normalization_206_gammaIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOp?assignvariableop_3_dense_model_103_batch_normalization_206_betaIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpFassignvariableop_4_dense_model_103_batch_normalization_206_moving_meanIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpJassignvariableop_5_dense_model_103_batch_normalization_206_moving_varianceIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp3assignvariableop_6_dense_model_103_dense_310_kernelIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp1assignvariableop_7_dense_model_103_dense_310_biasIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOp@assignvariableop_8_dense_model_103_batch_normalization_207_gammaIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp?assignvariableop_9_dense_model_103_batch_normalization_207_betaIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOpGassignvariableop_10_dense_model_103_batch_normalization_207_moving_meanIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpKassignvariableop_11_dense_model_103_batch_normalization_207_moving_varianceIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp4assignvariableop_12_dense_model_103_dense_311_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp2assignvariableop_13_dense_model_103_dense_311_biasIdentity_13:output:0"/device:CPU:0*&
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
AssignVariableOp_16AssignVariableOp;assignvariableop_16_adam_m_dense_model_103_dense_309_kernelIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOp;assignvariableop_17_adam_v_dense_model_103_dense_309_kernelIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp9assignvariableop_18_adam_m_dense_model_103_dense_309_biasIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp9assignvariableop_19_adam_v_dense_model_103_dense_309_biasIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOpHassignvariableop_20_adam_m_dense_model_103_batch_normalization_206_gammaIdentity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOpHassignvariableop_21_adam_v_dense_model_103_batch_normalization_206_gammaIdentity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOpGassignvariableop_22_adam_m_dense_model_103_batch_normalization_206_betaIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOpGassignvariableop_23_adam_v_dense_model_103_batch_normalization_206_betaIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOp;assignvariableop_24_adam_m_dense_model_103_dense_310_kernelIdentity_24:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOp;assignvariableop_25_adam_v_dense_model_103_dense_310_kernelIdentity_25:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOp9assignvariableop_26_adam_m_dense_model_103_dense_310_biasIdentity_26:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOp9assignvariableop_27_adam_v_dense_model_103_dense_310_biasIdentity_27:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpHassignvariableop_28_adam_m_dense_model_103_batch_normalization_207_gammaIdentity_28:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpHassignvariableop_29_adam_v_dense_model_103_batch_normalization_207_gammaIdentity_29:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOpGassignvariableop_30_adam_m_dense_model_103_batch_normalization_207_betaIdentity_30:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOpGassignvariableop_31_adam_v_dense_model_103_batch_normalization_207_betaIdentity_31:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOp;assignvariableop_32_adam_m_dense_model_103_dense_311_kernelIdentity_32:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOp;assignvariableop_33_adam_v_dense_model_103_dense_311_kernelIdentity_33:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOp9assignvariableop_34_adam_m_dense_model_103_dense_311_biasIdentity_34:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOp9assignvariableop_35_adam_v_dense_model_103_dense_311_biasIdentity_35:output:0"/device:CPU:0*&
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
_user_specified_name	total_1:E$A
?
_user_specified_name'%Adam/v/dense_model_103/dense_311/bias:E#A
?
_user_specified_name'%Adam/m/dense_model_103/dense_311/bias:G"C
A
_user_specified_name)'Adam/v/dense_model_103/dense_311/kernel:G!C
A
_user_specified_name)'Adam/m/dense_model_103/dense_311/kernel:S O
M
_user_specified_name53Adam/v/dense_model_103/batch_normalization_207/beta:SO
M
_user_specified_name53Adam/m/dense_model_103/batch_normalization_207/beta:TP
N
_user_specified_name64Adam/v/dense_model_103/batch_normalization_207/gamma:TP
N
_user_specified_name64Adam/m/dense_model_103/batch_normalization_207/gamma:EA
?
_user_specified_name'%Adam/v/dense_model_103/dense_310/bias:EA
?
_user_specified_name'%Adam/m/dense_model_103/dense_310/bias:GC
A
_user_specified_name)'Adam/v/dense_model_103/dense_310/kernel:GC
A
_user_specified_name)'Adam/m/dense_model_103/dense_310/kernel:SO
M
_user_specified_name53Adam/v/dense_model_103/batch_normalization_206/beta:SO
M
_user_specified_name53Adam/m/dense_model_103/batch_normalization_206/beta:TP
N
_user_specified_name64Adam/v/dense_model_103/batch_normalization_206/gamma:TP
N
_user_specified_name64Adam/m/dense_model_103/batch_normalization_206/gamma:EA
?
_user_specified_name'%Adam/v/dense_model_103/dense_309/bias:EA
?
_user_specified_name'%Adam/m/dense_model_103/dense_309/bias:GC
A
_user_specified_name)'Adam/v/dense_model_103/dense_309/kernel:GC
A
_user_specified_name)'Adam/m/dense_model_103/dense_309/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:>:
8
_user_specified_name dense_model_103/dense_311/bias:@<
:
_user_specified_name" dense_model_103/dense_311/kernel:WS
Q
_user_specified_name97dense_model_103/batch_normalization_207/moving_variance:SO
M
_user_specified_name53dense_model_103/batch_normalization_207/moving_mean:L
H
F
_user_specified_name.,dense_model_103/batch_normalization_207/beta:M	I
G
_user_specified_name/-dense_model_103/batch_normalization_207/gamma:>:
8
_user_specified_name dense_model_103/dense_310/bias:@<
:
_user_specified_name" dense_model_103/dense_310/kernel:WS
Q
_user_specified_name97dense_model_103/batch_normalization_206/moving_variance:SO
M
_user_specified_name53dense_model_103/batch_normalization_206/moving_mean:LH
F
_user_specified_name.,dense_model_103/batch_normalization_206/beta:MI
G
_user_specified_name/-dense_model_103/batch_normalization_206/gamma:>:
8
_user_specified_name dense_model_103/dense_309/bias:@<
:
_user_specified_name" dense_model_103/dense_309/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�%
�
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390489
input_1#
dense_309_390434:	�
dense_309_390436:,
batch_normalization_206_390439:,
batch_normalization_206_390441:,
batch_normalization_206_390443:,
batch_normalization_206_390445:"
dense_310_390458:
dense_310_390460:,
batch_normalization_207_390463:,
batch_normalization_207_390465:,
batch_normalization_207_390467:,
batch_normalization_207_390469:"
dense_311_390483:
dense_311_390485:
identity��/batch_normalization_206/StatefulPartitionedCall�/batch_normalization_207/StatefulPartitionedCall�!dense_309/StatefulPartitionedCall�!dense_310/StatefulPartitionedCall�!dense_311/StatefulPartitionedCall�
!dense_309/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_309_390434dense_309_390436*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_309_layer_call_and_return_conditional_losses_390433�
/batch_normalization_206/StatefulPartitionedCallStatefulPartitionedCall*dense_309/StatefulPartitionedCall:output:0batch_normalization_206_390439batch_normalization_206_390441batch_normalization_206_390443batch_normalization_206_390445*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390295�
!dense_310/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_206/StatefulPartitionedCall:output:0dense_310_390458dense_310_390460*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_310_layer_call_and_return_conditional_losses_390457�
/batch_normalization_207/StatefulPartitionedCallStatefulPartitionedCall*dense_310/StatefulPartitionedCall:output:0batch_normalization_207_390463batch_normalization_207_390465batch_normalization_207_390467batch_normalization_207_390469*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390375�
!dense_311/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_207/StatefulPartitionedCall:output:0dense_311_390483dense_311_390485*
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
E__inference_dense_311_layer_call_and_return_conditional_losses_390482y
IdentityIdentity*dense_311/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp0^batch_normalization_206/StatefulPartitionedCall0^batch_normalization_207/StatefulPartitionedCall"^dense_309/StatefulPartitionedCall"^dense_310/StatefulPartitionedCall"^dense_311/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2b
/batch_normalization_206/StatefulPartitionedCall/batch_normalization_206/StatefulPartitionedCall2b
/batch_normalization_207/StatefulPartitionedCall/batch_normalization_207/StatefulPartitionedCall2F
!dense_309/StatefulPartitionedCall!dense_309/StatefulPartitionedCall2F
!dense_310/StatefulPartitionedCall!dense_310/StatefulPartitionedCall2F
!dense_311/StatefulPartitionedCall!dense_311/StatefulPartitionedCall:&"
 
_user_specified_name390485:&"
 
_user_specified_name390483:&"
 
_user_specified_name390469:&"
 
_user_specified_name390467:&
"
 
_user_specified_name390465:&	"
 
_user_specified_name390463:&"
 
_user_specified_name390460:&"
 
_user_specified_name390458:&"
 
_user_specified_name390445:&"
 
_user_specified_name390443:&"
 
_user_specified_name390441:&"
 
_user_specified_name390439:&"
 
_user_specified_name390436:&"
 
_user_specified_name390434:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390395

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
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
:���������
 
_user_specified_nameinputs
�
�
$__inference_signature_wrapper_390655
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

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
!__inference__wrapped_model_390261o
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
 
_user_specified_name390651:&"
 
_user_specified_name390649:&"
 
_user_specified_name390647:&"
 
_user_specified_name390645:&
"
 
_user_specified_name390643:&	"
 
_user_specified_name390641:&"
 
_user_specified_name390639:&"
 
_user_specified_name390637:&"
 
_user_specified_name390635:&"
 
_user_specified_name390633:&"
 
_user_specified_name390631:&"
 
_user_specified_name390629:&"
 
_user_specified_name390627:&"
 
_user_specified_name390625:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1"�L
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
StatefulPartitionedCall:0���������tensorflow/serving/predict:��
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
0__inference_dense_model_103_layer_call_fn_390559
0__inference_dense_model_103_layer_call_fn_390592�
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
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390489
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390526�
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
!__inference__wrapped_model_390261input_1"�
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
3:1	�2 dense_model_103/dense_309/kernel
,:*2dense_model_103/dense_309/bias
;:92-dense_model_103/batch_normalization_206/gamma
::82,dense_model_103/batch_normalization_206/beta
C:A (23dense_model_103/batch_normalization_206/moving_mean
G:E (27dense_model_103/batch_normalization_206/moving_variance
2:02 dense_model_103/dense_310/kernel
,:*2dense_model_103/dense_310/bias
;:92-dense_model_103/batch_normalization_207/gamma
::82,dense_model_103/batch_normalization_207/beta
C:A (23dense_model_103/batch_normalization_207/moving_mean
G:E (27dense_model_103/batch_normalization_207/moving_variance
2:02 dense_model_103/dense_311/kernel
,:*2dense_model_103/dense_311/bias
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
0__inference_dense_model_103_layer_call_fn_390559input_1"�
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
0__inference_dense_model_103_layer_call_fn_390592input_1"�
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
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390489input_1"�
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
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390526input_1"�
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
*__inference_dense_309_layer_call_fn_390664�
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
E__inference_dense_309_layer_call_and_return_conditional_losses_390674�
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
8__inference_batch_normalization_206_layer_call_fn_390687
8__inference_batch_normalization_206_layer_call_fn_390700�
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
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390734
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390754�
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
*__inference_dense_310_layer_call_fn_390763�
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
E__inference_dense_310_layer_call_and_return_conditional_losses_390773�
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
8__inference_batch_normalization_207_layer_call_fn_390786
8__inference_batch_normalization_207_layer_call_fn_390799�
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
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390833
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390853�
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
*__inference_dense_311_layer_call_fn_390862�
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
E__inference_dense_311_layer_call_and_return_conditional_losses_390873�
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
$__inference_signature_wrapper_390655input_1"�
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
*__inference_dense_309_layer_call_fn_390664inputs"�
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
E__inference_dense_309_layer_call_and_return_conditional_losses_390674inputs"�
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
8__inference_batch_normalization_206_layer_call_fn_390687inputs"�
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
8__inference_batch_normalization_206_layer_call_fn_390700inputs"�
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
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390734inputs"�
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
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390754inputs"�
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
*__inference_dense_310_layer_call_fn_390763inputs"�
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
E__inference_dense_310_layer_call_and_return_conditional_losses_390773inputs"�
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
8__inference_batch_normalization_207_layer_call_fn_390786inputs"�
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
8__inference_batch_normalization_207_layer_call_fn_390799inputs"�
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
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390833inputs"�
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
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390853inputs"�
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
*__inference_dense_311_layer_call_fn_390862inputs"�
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
E__inference_dense_311_layer_call_and_return_conditional_losses_390873inputs"�
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
8:6	�2'Adam/m/dense_model_103/dense_309/kernel
8:6	�2'Adam/v/dense_model_103/dense_309/kernel
1:/2%Adam/m/dense_model_103/dense_309/bias
1:/2%Adam/v/dense_model_103/dense_309/bias
@:>24Adam/m/dense_model_103/batch_normalization_206/gamma
@:>24Adam/v/dense_model_103/batch_normalization_206/gamma
?:=23Adam/m/dense_model_103/batch_normalization_206/beta
?:=23Adam/v/dense_model_103/batch_normalization_206/beta
7:52'Adam/m/dense_model_103/dense_310/kernel
7:52'Adam/v/dense_model_103/dense_310/kernel
1:/2%Adam/m/dense_model_103/dense_310/bias
1:/2%Adam/v/dense_model_103/dense_310/bias
@:>24Adam/m/dense_model_103/batch_normalization_207/gamma
@:>24Adam/v/dense_model_103/batch_normalization_207/gamma
?:=23Adam/m/dense_model_103/batch_normalization_207/beta
?:=23Adam/v/dense_model_103/batch_normalization_207/beta
7:52'Adam/m/dense_model_103/dense_311/kernel
7:52'Adam/v/dense_model_103/dense_311/kernel
1:/2%Adam/m/dense_model_103/dense_311/bias
1:/2%Adam/v/dense_model_103/dense_311/bias
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
!__inference__wrapped_model_390261x1�.
'�$
"�
input_1����������
� "3�0
.
output_1"�
output_1����������
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390734m7�4
-�*
 �
inputs���������
p

 
� ",�)
"�
tensor_0���������
� �
S__inference_batch_normalization_206_layer_call_and_return_conditional_losses_390754m7�4
-�*
 �
inputs���������
p 

 
� ",�)
"�
tensor_0���������
� �
8__inference_batch_normalization_206_layer_call_fn_390687b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_206_layer_call_fn_390700b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390833m7�4
-�*
 �
inputs���������
p

 
� ",�)
"�
tensor_0���������
� �
S__inference_batch_normalization_207_layer_call_and_return_conditional_losses_390853m7�4
-�*
 �
inputs���������
p 

 
� ",�)
"�
tensor_0���������
� �
8__inference_batch_normalization_207_layer_call_fn_390786b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_207_layer_call_fn_390799b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
E__inference_dense_309_layer_call_and_return_conditional_losses_390674d0�-
&�#
!�
inputs����������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_309_layer_call_fn_390664Y0�-
&�#
!�
inputs����������
� "!�
unknown����������
E__inference_dense_310_layer_call_and_return_conditional_losses_390773c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_310_layer_call_fn_390763X/�,
%�"
 �
inputs���������
� "!�
unknown����������
E__inference_dense_311_layer_call_and_return_conditional_losses_390873c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_311_layer_call_fn_390862X/�,
%�"
 �
inputs���������
� "!�
unknown����������
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390489�A�>
'�$
"�
input_1����������
�

trainingp",�)
"�
tensor_0���������
� �
K__inference_dense_model_103_layer_call_and_return_conditional_losses_390526�A�>
'�$
"�
input_1����������
�

trainingp ",�)
"�
tensor_0���������
� �
0__inference_dense_model_103_layer_call_fn_390559vA�>
'�$
"�
input_1����������
�

trainingp"!�
unknown����������
0__inference_dense_model_103_layer_call_fn_390592vA�>
'�$
"�
input_1����������
�

trainingp "!�
unknown����������
$__inference_signature_wrapper_390655�<�9
� 
2�/
-
input_1"�
input_1����������"3�0
.
output_1"�
output_1���������