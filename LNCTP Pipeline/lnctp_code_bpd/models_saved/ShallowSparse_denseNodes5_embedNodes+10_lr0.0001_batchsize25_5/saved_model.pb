��
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
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
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
N

StringJoin
inputs*N

output"
Nint(0"
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
 �"serve*2.10.02v2.10.0-rc3-6-g359c3cdfc5f8��
�
1Adam/shallow_sparse_model_5124/dense_48049/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31Adam/shallow_sparse_model_5124/dense_48049/bias/v
�
EAdam/shallow_sparse_model_5124/dense_48049/bias/v/Read/ReadVariableOpReadVariableOp1Adam/shallow_sparse_model_5124/dense_48049/bias/v*
_output_shapes
:*
dtype0
�
3Adam/shallow_sparse_model_5124/dense_48049/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:Z*D
shared_name53Adam/shallow_sparse_model_5124/dense_48049/kernel/v
�
GAdam/shallow_sparse_model_5124/dense_48049/kernel/v/Read/ReadVariableOpReadVariableOp3Adam/shallow_sparse_model_5124/dense_48049/kernel/v*
_output_shapes

:Z*
dtype0
�
?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*P
shared_nameA?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/v
�
SAdam/shallow_sparse_model_5124/batch_normalization_32124/beta/v/Read/ReadVariableOpReadVariableOp?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/v*
_output_shapes
:Z*
dtype0
�
@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*Q
shared_nameB@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/v
�
TAdam/shallow_sparse_model_5124/batch_normalization_32124/gamma/v/Read/ReadVariableOpReadVariableOp@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/v*
_output_shapes
:Z*
dtype0
�
1Adam/shallow_sparse_model_5124/dense_48048/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*B
shared_name31Adam/shallow_sparse_model_5124/dense_48048/bias/v
�
EAdam/shallow_sparse_model_5124/dense_48048/bias/v/Read/ReadVariableOpReadVariableOp1Adam/shallow_sparse_model_5124/dense_48048/bias/v*
_output_shapes
:Z*
dtype0
�
3Adam/shallow_sparse_model_5124/dense_48048/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�Z*D
shared_name53Adam/shallow_sparse_model_5124/dense_48048/kernel/v
�
GAdam/shallow_sparse_model_5124/dense_48048/kernel/v/Read/ReadVariableOpReadVariableOp3Adam/shallow_sparse_model_5124/dense_48048/kernel/v*
_output_shapes
:	�Z*
dtype0
�
1Adam/shallow_sparse_model_5124/dense_48049/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*B
shared_name31Adam/shallow_sparse_model_5124/dense_48049/bias/m
�
EAdam/shallow_sparse_model_5124/dense_48049/bias/m/Read/ReadVariableOpReadVariableOp1Adam/shallow_sparse_model_5124/dense_48049/bias/m*
_output_shapes
:*
dtype0
�
3Adam/shallow_sparse_model_5124/dense_48049/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:Z*D
shared_name53Adam/shallow_sparse_model_5124/dense_48049/kernel/m
�
GAdam/shallow_sparse_model_5124/dense_48049/kernel/m/Read/ReadVariableOpReadVariableOp3Adam/shallow_sparse_model_5124/dense_48049/kernel/m*
_output_shapes

:Z*
dtype0
�
?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*P
shared_nameA?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/m
�
SAdam/shallow_sparse_model_5124/batch_normalization_32124/beta/m/Read/ReadVariableOpReadVariableOp?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/m*
_output_shapes
:Z*
dtype0
�
@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*Q
shared_nameB@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/m
�
TAdam/shallow_sparse_model_5124/batch_normalization_32124/gamma/m/Read/ReadVariableOpReadVariableOp@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/m*
_output_shapes
:Z*
dtype0
�
1Adam/shallow_sparse_model_5124/dense_48048/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*B
shared_name31Adam/shallow_sparse_model_5124/dense_48048/bias/m
�
EAdam/shallow_sparse_model_5124/dense_48048/bias/m/Read/ReadVariableOpReadVariableOp1Adam/shallow_sparse_model_5124/dense_48048/bias/m*
_output_shapes
:Z*
dtype0
�
3Adam/shallow_sparse_model_5124/dense_48048/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�Z*D
shared_name53Adam/shallow_sparse_model_5124/dense_48048/kernel/m
�
GAdam/shallow_sparse_model_5124/dense_48048/kernel/m/Read/ReadVariableOpReadVariableOp3Adam/shallow_sparse_model_5124/dense_48048/kernel/m*
_output_shapes
:	�Z*
dtype0
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
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
�
*shallow_sparse_model_5124/dense_48049/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*;
shared_name,*shallow_sparse_model_5124/dense_48049/bias
�
>shallow_sparse_model_5124/dense_48049/bias/Read/ReadVariableOpReadVariableOp*shallow_sparse_model_5124/dense_48049/bias*
_output_shapes
:*
dtype0
�
,shallow_sparse_model_5124/dense_48049/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:Z*=
shared_name.,shallow_sparse_model_5124/dense_48049/kernel
�
@shallow_sparse_model_5124/dense_48049/kernel/Read/ReadVariableOpReadVariableOp,shallow_sparse_model_5124/dense_48049/kernel*
_output_shapes

:Z*
dtype0
�
Cshallow_sparse_model_5124/batch_normalization_32124/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*T
shared_nameECshallow_sparse_model_5124/batch_normalization_32124/moving_variance
�
Wshallow_sparse_model_5124/batch_normalization_32124/moving_variance/Read/ReadVariableOpReadVariableOpCshallow_sparse_model_5124/batch_normalization_32124/moving_variance*
_output_shapes
:Z*
dtype0
�
?shallow_sparse_model_5124/batch_normalization_32124/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*P
shared_nameA?shallow_sparse_model_5124/batch_normalization_32124/moving_mean
�
Sshallow_sparse_model_5124/batch_normalization_32124/moving_mean/Read/ReadVariableOpReadVariableOp?shallow_sparse_model_5124/batch_normalization_32124/moving_mean*
_output_shapes
:Z*
dtype0
�
8shallow_sparse_model_5124/batch_normalization_32124/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*I
shared_name:8shallow_sparse_model_5124/batch_normalization_32124/beta
�
Lshallow_sparse_model_5124/batch_normalization_32124/beta/Read/ReadVariableOpReadVariableOp8shallow_sparse_model_5124/batch_normalization_32124/beta*
_output_shapes
:Z*
dtype0
�
9shallow_sparse_model_5124/batch_normalization_32124/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*J
shared_name;9shallow_sparse_model_5124/batch_normalization_32124/gamma
�
Mshallow_sparse_model_5124/batch_normalization_32124/gamma/Read/ReadVariableOpReadVariableOp9shallow_sparse_model_5124/batch_normalization_32124/gamma*
_output_shapes
:Z*
dtype0
�
*shallow_sparse_model_5124/dense_48048/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*;
shared_name,*shallow_sparse_model_5124/dense_48048/bias
�
>shallow_sparse_model_5124/dense_48048/bias/Read/ReadVariableOpReadVariableOp*shallow_sparse_model_5124/dense_48048/bias*
_output_shapes
:Z*
dtype0
�
,shallow_sparse_model_5124/dense_48048/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�Z*=
shared_name.,shallow_sparse_model_5124/dense_48048/kernel
�
@shallow_sparse_model_5124/dense_48048/kernel/Read/ReadVariableOpReadVariableOp,shallow_sparse_model_5124/dense_48048/kernel*
_output_shapes
:	�Z*
dtype0
|
serving_default_input_1Placeholder*(
_output_shapes
:����������*
dtype0*
shape:����������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1,shallow_sparse_model_5124/dense_48048/kernel*shallow_sparse_model_5124/dense_48048/bias?shallow_sparse_model_5124/batch_normalization_32124/moving_meanCshallow_sparse_model_5124/batch_normalization_32124/moving_variance8shallow_sparse_model_5124/batch_normalization_32124/beta9shallow_sparse_model_5124/batch_normalization_32124/gamma,shallow_sparse_model_5124/dense_48049/kernel*shallow_sparse_model_5124/dense_48049/bias*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� */
f*R(
&__inference_signature_wrapper_57495163

NoOpNoOp
�4
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�3
value�3B�3 B�3
�
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	embed
		norm1

out
	optimizer

signatures*
<
0
1
2
3
4
5
6
7*
.
0
1
2
3
4
5*
* 
�
non_trainable_variables

layers
metrics
layer_regularization_losses
layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
6
trace_0
trace_1
trace_2
trace_3* 
6
trace_0
trace_1
 trace_2
!trace_3* 
* 
�
"	variables
#trainable_variables
$regularization_losses
%	keras_api
&__call__
*'&call_and_return_all_conditional_losses

kernel
bias*
�
(	variables
)trainable_variables
*regularization_losses
+	keras_api
,__call__
*-&call_and_return_all_conditional_losses
.axis
	gamma
beta
moving_mean
moving_variance*
�
/	variables
0trainable_variables
1regularization_losses
2	keras_api
3__call__
*4&call_and_return_all_conditional_losses

kernel
bias*
�
5iter

6beta_1

7beta_2
	8decay
9learning_ratem]m^m_m`mambvcvdvevfvgvh*

:serving_default* 
lf
VARIABLE_VALUE,shallow_sparse_model_5124/dense_48048/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
jd
VARIABLE_VALUE*shallow_sparse_model_5124/dense_48048/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
ys
VARIABLE_VALUE9shallow_sparse_model_5124/batch_normalization_32124/gamma&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
xr
VARIABLE_VALUE8shallow_sparse_model_5124/batch_normalization_32124/beta&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE?shallow_sparse_model_5124/batch_normalization_32124/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
�}
VARIABLE_VALUECshallow_sparse_model_5124/batch_normalization_32124/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE,shallow_sparse_model_5124/dense_48049/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE*
jd
VARIABLE_VALUE*shallow_sparse_model_5124/dense_48049/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE*

0
1*

0
	1

2*

;0
<1*
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
0
1*

0
1*
* 
�
=non_trainable_variables

>layers
?metrics
@layer_regularization_losses
Alayer_metrics
"	variables
#trainable_variables
$regularization_losses
&__call__
*'&call_and_return_all_conditional_losses
&'"call_and_return_conditional_losses*

Btrace_0* 

Ctrace_0* 
 
0
1
2
3*

0
1*
* 
�
Dnon_trainable_variables

Elayers
Fmetrics
Glayer_regularization_losses
Hlayer_metrics
(	variables
)trainable_variables
*regularization_losses
,__call__
*-&call_and_return_all_conditional_losses
&-"call_and_return_conditional_losses*

Itrace_0
Jtrace_1* 

Ktrace_0
Ltrace_1* 
* 

0
1*

0
1*
* 
�
Mnon_trainable_variables

Nlayers
Ometrics
Player_regularization_losses
Qlayer_metrics
/	variables
0trainable_variables
1regularization_losses
3__call__
*4&call_and_return_all_conditional_losses
&4"call_and_return_conditional_losses*

Rtrace_0* 

Strace_0* 
LF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE*
PJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE*
PJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE*
NH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 
8
T	variables
U	keras_api
	Vtotal
	Wcount*
H
X	variables
Y	keras_api
	Ztotal
	[count
\
_fn_kwargs*
* 
* 
* 
* 
* 
* 
* 

0
1*
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
V0
W1*

T	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

Z0
[1*

X	variables*
SM
VARIABLE_VALUEtotal4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE*
* 
��
VARIABLE_VALUE3Adam/shallow_sparse_model_5124/dense_48048/kernel/mBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE1Adam/shallow_sparse_model_5124/dense_48048/bias/mBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/mBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/mBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE3Adam/shallow_sparse_model_5124/dense_48049/kernel/mBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE1Adam/shallow_sparse_model_5124/dense_48049/bias/mBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE3Adam/shallow_sparse_model_5124/dense_48048/kernel/vBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE1Adam/shallow_sparse_model_5124/dense_48048/bias/vBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/vBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/vBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE3Adam/shallow_sparse_model_5124/dense_48049/kernel/vBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE1Adam/shallow_sparse_model_5124/dense_48049/bias/vBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename@shallow_sparse_model_5124/dense_48048/kernel/Read/ReadVariableOp>shallow_sparse_model_5124/dense_48048/bias/Read/ReadVariableOpMshallow_sparse_model_5124/batch_normalization_32124/gamma/Read/ReadVariableOpLshallow_sparse_model_5124/batch_normalization_32124/beta/Read/ReadVariableOpSshallow_sparse_model_5124/batch_normalization_32124/moving_mean/Read/ReadVariableOpWshallow_sparse_model_5124/batch_normalization_32124/moving_variance/Read/ReadVariableOp@shallow_sparse_model_5124/dense_48049/kernel/Read/ReadVariableOp>shallow_sparse_model_5124/dense_48049/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpGAdam/shallow_sparse_model_5124/dense_48048/kernel/m/Read/ReadVariableOpEAdam/shallow_sparse_model_5124/dense_48048/bias/m/Read/ReadVariableOpTAdam/shallow_sparse_model_5124/batch_normalization_32124/gamma/m/Read/ReadVariableOpSAdam/shallow_sparse_model_5124/batch_normalization_32124/beta/m/Read/ReadVariableOpGAdam/shallow_sparse_model_5124/dense_48049/kernel/m/Read/ReadVariableOpEAdam/shallow_sparse_model_5124/dense_48049/bias/m/Read/ReadVariableOpGAdam/shallow_sparse_model_5124/dense_48048/kernel/v/Read/ReadVariableOpEAdam/shallow_sparse_model_5124/dense_48048/bias/v/Read/ReadVariableOpTAdam/shallow_sparse_model_5124/batch_normalization_32124/gamma/v/Read/ReadVariableOpSAdam/shallow_sparse_model_5124/batch_normalization_32124/beta/v/Read/ReadVariableOpGAdam/shallow_sparse_model_5124/dense_48049/kernel/v/Read/ReadVariableOpEAdam/shallow_sparse_model_5124/dense_48049/bias/v/Read/ReadVariableOpConst**
Tin#
!2	*
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
GPU 2J 8� **
f%R#
!__inference__traced_save_57495514
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename,shallow_sparse_model_5124/dense_48048/kernel*shallow_sparse_model_5124/dense_48048/bias9shallow_sparse_model_5124/batch_normalization_32124/gamma8shallow_sparse_model_5124/batch_normalization_32124/beta?shallow_sparse_model_5124/batch_normalization_32124/moving_meanCshallow_sparse_model_5124/batch_normalization_32124/moving_variance,shallow_sparse_model_5124/dense_48049/kernel*shallow_sparse_model_5124/dense_48049/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_ratetotal_1count_1totalcount3Adam/shallow_sparse_model_5124/dense_48048/kernel/m1Adam/shallow_sparse_model_5124/dense_48048/bias/m@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/m?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/m3Adam/shallow_sparse_model_5124/dense_48049/kernel/m1Adam/shallow_sparse_model_5124/dense_48049/bias/m3Adam/shallow_sparse_model_5124/dense_48048/kernel/v1Adam/shallow_sparse_model_5124/dense_48048/bias/v@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/v?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/v3Adam/shallow_sparse_model_5124/dense_48049/kernel/v1Adam/shallow_sparse_model_5124/dense_48049/bias/v*)
Tin"
 2*
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
GPU 2J 8� *-
f(R&
$__inference__traced_restore_57495611�
�	
�
<__inference_shallow_sparse_model_5124_layer_call_fn_57494982
input_1
unknown:	�Z
	unknown_0:Z
	unknown_1:Z
	unknown_2:Z
	unknown_3:Z
	unknown_4:Z
	unknown_5:Z
	unknown_6:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57494963o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
.__inference_dense_48049_layer_call_fn_57495393

inputs
unknown:Z
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
GPU 2J 8� *R
fMRK
I__inference_dense_48049_layer_call_and_return_conditional_losses_57494956o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������Z: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�
�
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495350

inputs*
cast_readvariableop_resource:Z,
cast_1_readvariableop_resource:Z,
cast_2_readvariableop_resource:Z,
cast_3_readvariableop_resource:Z
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:Z*
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
:ZP
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:Zm
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:Zc
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������Zk
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:Zm
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:Zr
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Zb
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������Z�
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������Z: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�$
�
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494902

inputs5
'assignmovingavg_readvariableop_resource:Z7
)assignmovingavg_1_readvariableop_resource:Z*
cast_readvariableop_resource:Z,
cast_1_readvariableop_resource:Z
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:Z�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������Zl
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:Z*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:Z*
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
:Z*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:Zx
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:Z�
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
:Z*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:Z~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:Z�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:Z*
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
:ZP
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:Zm
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:Zc
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������Zh
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:Zk
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:Zr
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Zb
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������Z�
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������Z: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�	
�
I__inference_dense_48048_layer_call_and_return_conditional_losses_57495304

inputs1
matmul_readvariableop_resource:	�Z-
biasadd_readvariableop_resource:Z
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Zr
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������Zw
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57494963
x'
dense_48048_57494931:	�Z"
dense_48048_57494933:Z0
"batch_normalization_32124_57494936:Z0
"batch_normalization_32124_57494938:Z0
"batch_normalization_32124_57494940:Z0
"batch_normalization_32124_57494942:Z&
dense_48049_57494957:Z"
dense_48049_57494959:
identity��1batch_normalization_32124/StatefulPartitionedCall�#dense_48048/StatefulPartitionedCall�#dense_48049/StatefulPartitionedCall�
#dense_48048/StatefulPartitionedCallStatefulPartitionedCallxdense_48048_57494931dense_48048_57494933*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_48048_layer_call_and_return_conditional_losses_57494930�
1batch_normalization_32124/StatefulPartitionedCallStatefulPartitionedCall,dense_48048/StatefulPartitionedCall:output:0"batch_normalization_32124_57494936"batch_normalization_32124_57494938"batch_normalization_32124_57494940"batch_normalization_32124_57494942*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494855�
#dense_48049/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_32124/StatefulPartitionedCall:output:0dense_48049_57494957dense_48049_57494959*
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
GPU 2J 8� *R
fMRK
I__inference_dense_48049_layer_call_and_return_conditional_losses_57494956{
IdentityIdentity,dense_48049/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_32124/StatefulPartitionedCall$^dense_48048/StatefulPartitionedCall$^dense_48049/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_32124/StatefulPartitionedCall1batch_normalization_32124/StatefulPartitionedCall2J
#dense_48048/StatefulPartitionedCall#dense_48048/StatefulPartitionedCall2J
#dense_48049/StatefulPartitionedCall#dense_48049/StatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�
�
.__inference_dense_48048_layer_call_fn_57495294

inputs
unknown:	�Z
	unknown_0:Z
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_48048_layer_call_and_return_conditional_losses_57494930o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������Z`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495048
x'
dense_48048_57495028:	�Z"
dense_48048_57495030:Z0
"batch_normalization_32124_57495033:Z0
"batch_normalization_32124_57495035:Z0
"batch_normalization_32124_57495037:Z0
"batch_normalization_32124_57495039:Z&
dense_48049_57495042:Z"
dense_48049_57495044:
identity��1batch_normalization_32124/StatefulPartitionedCall�#dense_48048/StatefulPartitionedCall�#dense_48049/StatefulPartitionedCall�
#dense_48048/StatefulPartitionedCallStatefulPartitionedCallxdense_48048_57495028dense_48048_57495030*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_48048_layer_call_and_return_conditional_losses_57494930�
1batch_normalization_32124/StatefulPartitionedCallStatefulPartitionedCall,dense_48048/StatefulPartitionedCall:output:0"batch_normalization_32124_57495033"batch_normalization_32124_57495035"batch_normalization_32124_57495037"batch_normalization_32124_57495039*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494902�
#dense_48049/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_32124/StatefulPartitionedCall:output:0dense_48049_57495042dense_48049_57495044*
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
GPU 2J 8� *R
fMRK
I__inference_dense_48049_layer_call_and_return_conditional_losses_57494956{
IdentityIdentity,dense_48049/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_32124/StatefulPartitionedCall$^dense_48048/StatefulPartitionedCall$^dense_48049/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_32124/StatefulPartitionedCall1batch_normalization_32124/StatefulPartitionedCall2J
#dense_48048/StatefulPartitionedCall#dense_48048/StatefulPartitionedCall2J
#dense_48049/StatefulPartitionedCall#dense_48049/StatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�	
�
<__inference_shallow_sparse_model_5124_layer_call_fn_57495088
input_1
unknown:	�Z
	unknown_0:Z
	unknown_1:Z
	unknown_2:Z
	unknown_3:Z
	unknown_4:Z
	unknown_5:Z
	unknown_6:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495048o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�	
�
<__inference_shallow_sparse_model_5124_layer_call_fn_57495205
x
unknown:	�Z
	unknown_0:Z
	unknown_1:Z
	unknown_2:Z
	unknown_3:Z
	unknown_4:Z
	unknown_5:Z
	unknown_6:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495048o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�	
�
&__inference_signature_wrapper_57495163
input_1
unknown:	�Z
	unknown_0:Z
	unknown_1:Z
	unknown_2:Z
	unknown_3:Z
	unknown_4:Z
	unknown_5:Z
	unknown_6:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *,
f'R%
#__inference__wrapped_model_57494831o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
<__inference_batch_normalization_32124_layer_call_fn_57495317

inputs
unknown:Z
	unknown_0:Z
	unknown_1:Z
	unknown_2:Z
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494855o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������Z`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������Z: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�?
�

#__inference__wrapped_model_57494831
input_1W
Dshallow_sparse_model_5124_dense_48048_matmul_readvariableop_resource:	�ZS
Eshallow_sparse_model_5124_dense_48048_biasadd_readvariableop_resource:Z^
Pshallow_sparse_model_5124_batch_normalization_32124_cast_readvariableop_resource:Z`
Rshallow_sparse_model_5124_batch_normalization_32124_cast_1_readvariableop_resource:Z`
Rshallow_sparse_model_5124_batch_normalization_32124_cast_2_readvariableop_resource:Z`
Rshallow_sparse_model_5124_batch_normalization_32124_cast_3_readvariableop_resource:ZV
Dshallow_sparse_model_5124_dense_48049_matmul_readvariableop_resource:ZS
Eshallow_sparse_model_5124_dense_48049_biasadd_readvariableop_resource:
identity��Gshallow_sparse_model_5124/batch_normalization_32124/Cast/ReadVariableOp�Ishallow_sparse_model_5124/batch_normalization_32124/Cast_1/ReadVariableOp�Ishallow_sparse_model_5124/batch_normalization_32124/Cast_2/ReadVariableOp�Ishallow_sparse_model_5124/batch_normalization_32124/Cast_3/ReadVariableOp�<shallow_sparse_model_5124/dense_48048/BiasAdd/ReadVariableOp�;shallow_sparse_model_5124/dense_48048/MatMul/ReadVariableOp�<shallow_sparse_model_5124/dense_48049/BiasAdd/ReadVariableOp�;shallow_sparse_model_5124/dense_48049/MatMul/ReadVariableOp�
;shallow_sparse_model_5124/dense_48048/MatMul/ReadVariableOpReadVariableOpDshallow_sparse_model_5124_dense_48048_matmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0�
,shallow_sparse_model_5124/dense_48048/MatMulMatMulinput_1Cshallow_sparse_model_5124/dense_48048/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
<shallow_sparse_model_5124/dense_48048/BiasAdd/ReadVariableOpReadVariableOpEshallow_sparse_model_5124_dense_48048_biasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0�
-shallow_sparse_model_5124/dense_48048/BiasAddBiasAdd6shallow_sparse_model_5124/dense_48048/MatMul:product:0Dshallow_sparse_model_5124/dense_48048/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
Gshallow_sparse_model_5124/batch_normalization_32124/Cast/ReadVariableOpReadVariableOpPshallow_sparse_model_5124_batch_normalization_32124_cast_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Ishallow_sparse_model_5124/batch_normalization_32124/Cast_1/ReadVariableOpReadVariableOpRshallow_sparse_model_5124_batch_normalization_32124_cast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Ishallow_sparse_model_5124/batch_normalization_32124/Cast_2/ReadVariableOpReadVariableOpRshallow_sparse_model_5124_batch_normalization_32124_cast_2_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Ishallow_sparse_model_5124/batch_normalization_32124/Cast_3/ReadVariableOpReadVariableOpRshallow_sparse_model_5124_batch_normalization_32124_cast_3_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Cshallow_sparse_model_5124/batch_normalization_32124/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
Ashallow_sparse_model_5124/batch_normalization_32124/batchnorm/addAddV2Qshallow_sparse_model_5124/batch_normalization_32124/Cast_1/ReadVariableOp:value:0Lshallow_sparse_model_5124/batch_normalization_32124/batchnorm/add/y:output:0*
T0*
_output_shapes
:Z�
Cshallow_sparse_model_5124/batch_normalization_32124/batchnorm/RsqrtRsqrtEshallow_sparse_model_5124/batch_normalization_32124/batchnorm/add:z:0*
T0*
_output_shapes
:Z�
Ashallow_sparse_model_5124/batch_normalization_32124/batchnorm/mulMulGshallow_sparse_model_5124/batch_normalization_32124/batchnorm/Rsqrt:y:0Qshallow_sparse_model_5124/batch_normalization_32124/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:Z�
Cshallow_sparse_model_5124/batch_normalization_32124/batchnorm/mul_1Mul6shallow_sparse_model_5124/dense_48048/BiasAdd:output:0Eshallow_sparse_model_5124/batch_normalization_32124/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������Z�
Cshallow_sparse_model_5124/batch_normalization_32124/batchnorm/mul_2MulOshallow_sparse_model_5124/batch_normalization_32124/Cast/ReadVariableOp:value:0Eshallow_sparse_model_5124/batch_normalization_32124/batchnorm/mul:z:0*
T0*
_output_shapes
:Z�
Ashallow_sparse_model_5124/batch_normalization_32124/batchnorm/subSubQshallow_sparse_model_5124/batch_normalization_32124/Cast_2/ReadVariableOp:value:0Gshallow_sparse_model_5124/batch_normalization_32124/batchnorm/mul_2:z:0*
T0*
_output_shapes
:Z�
Cshallow_sparse_model_5124/batch_normalization_32124/batchnorm/add_1AddV2Gshallow_sparse_model_5124/batch_normalization_32124/batchnorm/mul_1:z:0Eshallow_sparse_model_5124/batch_normalization_32124/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Z�
;shallow_sparse_model_5124/dense_48049/MatMul/ReadVariableOpReadVariableOpDshallow_sparse_model_5124_dense_48049_matmul_readvariableop_resource*
_output_shapes

:Z*
dtype0�
,shallow_sparse_model_5124/dense_48049/MatMulMatMulGshallow_sparse_model_5124/batch_normalization_32124/batchnorm/add_1:z:0Cshallow_sparse_model_5124/dense_48049/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
<shallow_sparse_model_5124/dense_48049/BiasAdd/ReadVariableOpReadVariableOpEshallow_sparse_model_5124_dense_48049_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
-shallow_sparse_model_5124/dense_48049/BiasAddBiasAdd6shallow_sparse_model_5124/dense_48049/MatMul:product:0Dshallow_sparse_model_5124/dense_48049/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
-shallow_sparse_model_5124/dense_48049/SigmoidSigmoid6shallow_sparse_model_5124/dense_48049/BiasAdd:output:0*
T0*'
_output_shapes
:����������
IdentityIdentity1shallow_sparse_model_5124/dense_48049/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOpH^shallow_sparse_model_5124/batch_normalization_32124/Cast/ReadVariableOpJ^shallow_sparse_model_5124/batch_normalization_32124/Cast_1/ReadVariableOpJ^shallow_sparse_model_5124/batch_normalization_32124/Cast_2/ReadVariableOpJ^shallow_sparse_model_5124/batch_normalization_32124/Cast_3/ReadVariableOp=^shallow_sparse_model_5124/dense_48048/BiasAdd/ReadVariableOp<^shallow_sparse_model_5124/dense_48048/MatMul/ReadVariableOp=^shallow_sparse_model_5124/dense_48049/BiasAdd/ReadVariableOp<^shallow_sparse_model_5124/dense_48049/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2�
Gshallow_sparse_model_5124/batch_normalization_32124/Cast/ReadVariableOpGshallow_sparse_model_5124/batch_normalization_32124/Cast/ReadVariableOp2�
Ishallow_sparse_model_5124/batch_normalization_32124/Cast_1/ReadVariableOpIshallow_sparse_model_5124/batch_normalization_32124/Cast_1/ReadVariableOp2�
Ishallow_sparse_model_5124/batch_normalization_32124/Cast_2/ReadVariableOpIshallow_sparse_model_5124/batch_normalization_32124/Cast_2/ReadVariableOp2�
Ishallow_sparse_model_5124/batch_normalization_32124/Cast_3/ReadVariableOpIshallow_sparse_model_5124/batch_normalization_32124/Cast_3/ReadVariableOp2|
<shallow_sparse_model_5124/dense_48048/BiasAdd/ReadVariableOp<shallow_sparse_model_5124/dense_48048/BiasAdd/ReadVariableOp2z
;shallow_sparse_model_5124/dense_48048/MatMul/ReadVariableOp;shallow_sparse_model_5124/dense_48048/MatMul/ReadVariableOp2|
<shallow_sparse_model_5124/dense_48049/BiasAdd/ReadVariableOp<shallow_sparse_model_5124/dense_48049/BiasAdd/ReadVariableOp2z
;shallow_sparse_model_5124/dense_48049/MatMul/ReadVariableOp;shallow_sparse_model_5124/dense_48049/MatMul/ReadVariableOp:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�

�
I__inference_dense_48049_layer_call_and_return_conditional_losses_57495404

inputs0
matmul_readvariableop_resource:Z-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:Z*
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
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������Z: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�|
�
$__inference__traced_restore_57495611
file_prefixP
=assignvariableop_shallow_sparse_model_5124_dense_48048_kernel:	�ZK
=assignvariableop_1_shallow_sparse_model_5124_dense_48048_bias:ZZ
Lassignvariableop_2_shallow_sparse_model_5124_batch_normalization_32124_gamma:ZY
Kassignvariableop_3_shallow_sparse_model_5124_batch_normalization_32124_beta:Z`
Rassignvariableop_4_shallow_sparse_model_5124_batch_normalization_32124_moving_mean:Zd
Vassignvariableop_5_shallow_sparse_model_5124_batch_normalization_32124_moving_variance:ZQ
?assignvariableop_6_shallow_sparse_model_5124_dense_48049_kernel:ZK
=assignvariableop_7_shallow_sparse_model_5124_dense_48049_bias:&
assignvariableop_8_adam_iter:	 (
assignvariableop_9_adam_beta_1: )
assignvariableop_10_adam_beta_2: (
assignvariableop_11_adam_decay: 0
&assignvariableop_12_adam_learning_rate: %
assignvariableop_13_total_1: %
assignvariableop_14_count_1: #
assignvariableop_15_total: #
assignvariableop_16_count: Z
Gassignvariableop_17_adam_shallow_sparse_model_5124_dense_48048_kernel_m:	�ZS
Eassignvariableop_18_adam_shallow_sparse_model_5124_dense_48048_bias_m:Zb
Tassignvariableop_19_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_m:Za
Sassignvariableop_20_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_m:ZY
Gassignvariableop_21_adam_shallow_sparse_model_5124_dense_48049_kernel_m:ZS
Eassignvariableop_22_adam_shallow_sparse_model_5124_dense_48049_bias_m:Z
Gassignvariableop_23_adam_shallow_sparse_model_5124_dense_48048_kernel_v:	�ZS
Eassignvariableop_24_adam_shallow_sparse_model_5124_dense_48048_bias_v:Zb
Tassignvariableop_25_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_v:Za
Sassignvariableop_26_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_v:ZY
Gassignvariableop_27_adam_shallow_sparse_model_5124_dense_48049_kernel_v:ZS
Eassignvariableop_28_adam_shallow_sparse_model_5124_dense_48049_bias_v:
identity_30��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_3�AssignVariableOp_4�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*�
value�B�B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*O
valueFBDB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapesz
x::::::::::::::::::::::::::::::*,
dtypes"
 2	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOp=assignvariableop_shallow_sparse_model_5124_dense_48048_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp=assignvariableop_1_shallow_sparse_model_5124_dense_48048_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOpLassignvariableop_2_shallow_sparse_model_5124_batch_normalization_32124_gammaIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpKassignvariableop_3_shallow_sparse_model_5124_batch_normalization_32124_betaIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpRassignvariableop_4_shallow_sparse_model_5124_batch_normalization_32124_moving_meanIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpVassignvariableop_5_shallow_sparse_model_5124_batch_normalization_32124_moving_varianceIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp?assignvariableop_6_shallow_sparse_model_5124_dense_48049_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp=assignvariableop_7_shallow_sparse_model_5124_dense_48049_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_8AssignVariableOpassignvariableop_8_adam_iterIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOpassignvariableop_9_adam_beta_1Identity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOpassignvariableop_10_adam_beta_2Identity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpassignvariableop_11_adam_decayIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp&assignvariableop_12_adam_learning_rateIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOpassignvariableop_13_total_1Identity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOpassignvariableop_14_count_1Identity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOpassignvariableop_15_totalIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOpassignvariableop_16_countIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOpGassignvariableop_17_adam_shallow_sparse_model_5124_dense_48048_kernel_mIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOpEassignvariableop_18_adam_shallow_sparse_model_5124_dense_48048_bias_mIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOpTassignvariableop_19_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_mIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOpSassignvariableop_20_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_mIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOpGassignvariableop_21_adam_shallow_sparse_model_5124_dense_48049_kernel_mIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOpEassignvariableop_22_adam_shallow_sparse_model_5124_dense_48049_bias_mIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOpGassignvariableop_23_adam_shallow_sparse_model_5124_dense_48048_kernel_vIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOpEassignvariableop_24_adam_shallow_sparse_model_5124_dense_48048_bias_vIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOpTassignvariableop_25_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_vIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpSassignvariableop_26_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_vIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOpGassignvariableop_27_adam_shallow_sparse_model_5124_dense_48049_kernel_vIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpEassignvariableop_28_adam_shallow_sparse_model_5124_dense_48049_bias_vIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype01
NoOpNoOp"/device:CPU:0*
_output_shapes
 �
Identity_29Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_30IdentityIdentity_29:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_30Identity_30:output:0*O
_input_shapes>
<: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
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
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�	
�
I__inference_dense_48048_layer_call_and_return_conditional_losses_57494930

inputs1
matmul_readvariableop_resource:	�Z-
biasadd_readvariableop_resource:Z
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Zr
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������Zw
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495111
input_1'
dense_48048_57495091:	�Z"
dense_48048_57495093:Z0
"batch_normalization_32124_57495096:Z0
"batch_normalization_32124_57495098:Z0
"batch_normalization_32124_57495100:Z0
"batch_normalization_32124_57495102:Z&
dense_48049_57495105:Z"
dense_48049_57495107:
identity��1batch_normalization_32124/StatefulPartitionedCall�#dense_48048/StatefulPartitionedCall�#dense_48049/StatefulPartitionedCall�
#dense_48048/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_48048_57495091dense_48048_57495093*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_48048_layer_call_and_return_conditional_losses_57494930�
1batch_normalization_32124/StatefulPartitionedCallStatefulPartitionedCall,dense_48048/StatefulPartitionedCall:output:0"batch_normalization_32124_57495096"batch_normalization_32124_57495098"batch_normalization_32124_57495100"batch_normalization_32124_57495102*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494855�
#dense_48049/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_32124/StatefulPartitionedCall:output:0dense_48049_57495105dense_48049_57495107*
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
GPU 2J 8� *R
fMRK
I__inference_dense_48049_layer_call_and_return_conditional_losses_57494956{
IdentityIdentity,dense_48049/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_32124/StatefulPartitionedCall$^dense_48048/StatefulPartitionedCall$^dense_48049/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_32124/StatefulPartitionedCall1batch_normalization_32124/StatefulPartitionedCall2J
#dense_48048/StatefulPartitionedCall#dense_48048/StatefulPartitionedCall2J
#dense_48049/StatefulPartitionedCall#dense_48049/StatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�H
�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495285
x=
*dense_48048_matmul_readvariableop_resource:	�Z9
+dense_48048_biasadd_readvariableop_resource:ZO
Abatch_normalization_32124_assignmovingavg_readvariableop_resource:ZQ
Cbatch_normalization_32124_assignmovingavg_1_readvariableop_resource:ZD
6batch_normalization_32124_cast_readvariableop_resource:ZF
8batch_normalization_32124_cast_1_readvariableop_resource:Z<
*dense_48049_matmul_readvariableop_resource:Z9
+dense_48049_biasadd_readvariableop_resource:
identity��)batch_normalization_32124/AssignMovingAvg�8batch_normalization_32124/AssignMovingAvg/ReadVariableOp�+batch_normalization_32124/AssignMovingAvg_1�:batch_normalization_32124/AssignMovingAvg_1/ReadVariableOp�-batch_normalization_32124/Cast/ReadVariableOp�/batch_normalization_32124/Cast_1/ReadVariableOp�"dense_48048/BiasAdd/ReadVariableOp�!dense_48048/MatMul/ReadVariableOp�"dense_48049/BiasAdd/ReadVariableOp�!dense_48049/MatMul/ReadVariableOp�
!dense_48048/MatMul/ReadVariableOpReadVariableOp*dense_48048_matmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0|
dense_48048/MatMulMatMulx)dense_48048/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
"dense_48048/BiasAdd/ReadVariableOpReadVariableOp+dense_48048_biasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0�
dense_48048/BiasAddBiasAdddense_48048/MatMul:product:0*dense_48048/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
8batch_normalization_32124/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
&batch_normalization_32124/moments/meanMeandense_48048/BiasAdd:output:0Abatch_normalization_32124/moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(�
.batch_normalization_32124/moments/StopGradientStopGradient/batch_normalization_32124/moments/mean:output:0*
T0*
_output_shapes

:Z�
3batch_normalization_32124/moments/SquaredDifferenceSquaredDifferencedense_48048/BiasAdd:output:07batch_normalization_32124/moments/StopGradient:output:0*
T0*'
_output_shapes
:���������Z�
<batch_normalization_32124/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
*batch_normalization_32124/moments/varianceMean7batch_normalization_32124/moments/SquaredDifference:z:0Ebatch_normalization_32124/moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(�
)batch_normalization_32124/moments/SqueezeSqueeze/batch_normalization_32124/moments/mean:output:0*
T0*
_output_shapes
:Z*
squeeze_dims
 �
+batch_normalization_32124/moments/Squeeze_1Squeeze3batch_normalization_32124/moments/variance:output:0*
T0*
_output_shapes
:Z*
squeeze_dims
 t
/batch_normalization_32124/AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
8batch_normalization_32124/AssignMovingAvg/ReadVariableOpReadVariableOpAbatch_normalization_32124_assignmovingavg_readvariableop_resource*
_output_shapes
:Z*
dtype0�
-batch_normalization_32124/AssignMovingAvg/subSub@batch_normalization_32124/AssignMovingAvg/ReadVariableOp:value:02batch_normalization_32124/moments/Squeeze:output:0*
T0*
_output_shapes
:Z�
-batch_normalization_32124/AssignMovingAvg/mulMul1batch_normalization_32124/AssignMovingAvg/sub:z:08batch_normalization_32124/AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:Z�
)batch_normalization_32124/AssignMovingAvgAssignSubVariableOpAbatch_normalization_32124_assignmovingavg_readvariableop_resource1batch_normalization_32124/AssignMovingAvg/mul:z:09^batch_normalization_32124/AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0v
1batch_normalization_32124/AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
:batch_normalization_32124/AssignMovingAvg_1/ReadVariableOpReadVariableOpCbatch_normalization_32124_assignmovingavg_1_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_32124/AssignMovingAvg_1/subSubBbatch_normalization_32124/AssignMovingAvg_1/ReadVariableOp:value:04batch_normalization_32124/moments/Squeeze_1:output:0*
T0*
_output_shapes
:Z�
/batch_normalization_32124/AssignMovingAvg_1/mulMul3batch_normalization_32124/AssignMovingAvg_1/sub:z:0:batch_normalization_32124/AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:Z�
+batch_normalization_32124/AssignMovingAvg_1AssignSubVariableOpCbatch_normalization_32124_assignmovingavg_1_readvariableop_resource3batch_normalization_32124/AssignMovingAvg_1/mul:z:0;^batch_normalization_32124/AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0�
-batch_normalization_32124/Cast/ReadVariableOpReadVariableOp6batch_normalization_32124_cast_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_32124/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_32124_cast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0n
)batch_normalization_32124/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_32124/batchnorm/addAddV24batch_normalization_32124/moments/Squeeze_1:output:02batch_normalization_32124/batchnorm/add/y:output:0*
T0*
_output_shapes
:Z�
)batch_normalization_32124/batchnorm/RsqrtRsqrt+batch_normalization_32124/batchnorm/add:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_32124/batchnorm/mulMul-batch_normalization_32124/batchnorm/Rsqrt:y:07batch_normalization_32124/Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:Z�
)batch_normalization_32124/batchnorm/mul_1Muldense_48048/BiasAdd:output:0+batch_normalization_32124/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������Z�
)batch_normalization_32124/batchnorm/mul_2Mul2batch_normalization_32124/moments/Squeeze:output:0+batch_normalization_32124/batchnorm/mul:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_32124/batchnorm/subSub5batch_normalization_32124/Cast/ReadVariableOp:value:0-batch_normalization_32124/batchnorm/mul_2:z:0*
T0*
_output_shapes
:Z�
)batch_normalization_32124/batchnorm/add_1AddV2-batch_normalization_32124/batchnorm/mul_1:z:0+batch_normalization_32124/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Z�
!dense_48049/MatMul/ReadVariableOpReadVariableOp*dense_48049_matmul_readvariableop_resource*
_output_shapes

:Z*
dtype0�
dense_48049/MatMulMatMul-batch_normalization_32124/batchnorm/add_1:z:0)dense_48049/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
"dense_48049/BiasAdd/ReadVariableOpReadVariableOp+dense_48049_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_48049/BiasAddBiasAdddense_48049/MatMul:product:0*dense_48049/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������n
dense_48049/SigmoidSigmoiddense_48049/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
IdentityIdentitydense_48049/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp*^batch_normalization_32124/AssignMovingAvg9^batch_normalization_32124/AssignMovingAvg/ReadVariableOp,^batch_normalization_32124/AssignMovingAvg_1;^batch_normalization_32124/AssignMovingAvg_1/ReadVariableOp.^batch_normalization_32124/Cast/ReadVariableOp0^batch_normalization_32124/Cast_1/ReadVariableOp#^dense_48048/BiasAdd/ReadVariableOp"^dense_48048/MatMul/ReadVariableOp#^dense_48049/BiasAdd/ReadVariableOp"^dense_48049/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2V
)batch_normalization_32124/AssignMovingAvg)batch_normalization_32124/AssignMovingAvg2t
8batch_normalization_32124/AssignMovingAvg/ReadVariableOp8batch_normalization_32124/AssignMovingAvg/ReadVariableOp2Z
+batch_normalization_32124/AssignMovingAvg_1+batch_normalization_32124/AssignMovingAvg_12x
:batch_normalization_32124/AssignMovingAvg_1/ReadVariableOp:batch_normalization_32124/AssignMovingAvg_1/ReadVariableOp2^
-batch_normalization_32124/Cast/ReadVariableOp-batch_normalization_32124/Cast/ReadVariableOp2b
/batch_normalization_32124/Cast_1/ReadVariableOp/batch_normalization_32124/Cast_1/ReadVariableOp2H
"dense_48048/BiasAdd/ReadVariableOp"dense_48048/BiasAdd/ReadVariableOp2F
!dense_48048/MatMul/ReadVariableOp!dense_48048/MatMul/ReadVariableOp2H
"dense_48049/BiasAdd/ReadVariableOp"dense_48049/BiasAdd/ReadVariableOp2F
!dense_48049/MatMul/ReadVariableOp!dense_48049/MatMul/ReadVariableOp:K G
(
_output_shapes
:����������

_user_specified_namex
�$
�
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495384

inputs5
'assignmovingavg_readvariableop_resource:Z7
)assignmovingavg_1_readvariableop_resource:Z*
cast_readvariableop_resource:Z,
cast_1_readvariableop_resource:Z
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:Z�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������Zl
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:Z*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:Z*
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
:Z*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:Zx
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:Z�
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
:Z*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:Z~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:Z�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:Z*
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
:ZP
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:Zm
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:Zc
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������Zh
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:Zk
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:Zr
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Zb
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������Z�
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������Z: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�G
�
!__inference__traced_save_57495514
file_prefixK
Gsavev2_shallow_sparse_model_5124_dense_48048_kernel_read_readvariableopI
Esavev2_shallow_sparse_model_5124_dense_48048_bias_read_readvariableopX
Tsavev2_shallow_sparse_model_5124_batch_normalization_32124_gamma_read_readvariableopW
Ssavev2_shallow_sparse_model_5124_batch_normalization_32124_beta_read_readvariableop^
Zsavev2_shallow_sparse_model_5124_batch_normalization_32124_moving_mean_read_readvariableopb
^savev2_shallow_sparse_model_5124_batch_normalization_32124_moving_variance_read_readvariableopK
Gsavev2_shallow_sparse_model_5124_dense_48049_kernel_read_readvariableopI
Esavev2_shallow_sparse_model_5124_dense_48049_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableopR
Nsavev2_adam_shallow_sparse_model_5124_dense_48048_kernel_m_read_readvariableopP
Lsavev2_adam_shallow_sparse_model_5124_dense_48048_bias_m_read_readvariableop_
[savev2_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_m_read_readvariableop^
Zsavev2_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_m_read_readvariableopR
Nsavev2_adam_shallow_sparse_model_5124_dense_48049_kernel_m_read_readvariableopP
Lsavev2_adam_shallow_sparse_model_5124_dense_48049_bias_m_read_readvariableopR
Nsavev2_adam_shallow_sparse_model_5124_dense_48048_kernel_v_read_readvariableopP
Lsavev2_adam_shallow_sparse_model_5124_dense_48048_bias_v_read_readvariableop_
[savev2_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_v_read_readvariableop^
Zsavev2_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_v_read_readvariableopR
Nsavev2_adam_shallow_sparse_model_5124_dense_48049_kernel_v_read_readvariableopP
Lsavev2_adam_shallow_sparse_model_5124_dense_48049_bias_v_read_readvariableop
savev2_const

identity_1��MergeV2Checkpointsw
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
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*�
value�B�B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*O
valueFBDB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0Gsavev2_shallow_sparse_model_5124_dense_48048_kernel_read_readvariableopEsavev2_shallow_sparse_model_5124_dense_48048_bias_read_readvariableopTsavev2_shallow_sparse_model_5124_batch_normalization_32124_gamma_read_readvariableopSsavev2_shallow_sparse_model_5124_batch_normalization_32124_beta_read_readvariableopZsavev2_shallow_sparse_model_5124_batch_normalization_32124_moving_mean_read_readvariableop^savev2_shallow_sparse_model_5124_batch_normalization_32124_moving_variance_read_readvariableopGsavev2_shallow_sparse_model_5124_dense_48049_kernel_read_readvariableopEsavev2_shallow_sparse_model_5124_dense_48049_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopNsavev2_adam_shallow_sparse_model_5124_dense_48048_kernel_m_read_readvariableopLsavev2_adam_shallow_sparse_model_5124_dense_48048_bias_m_read_readvariableop[savev2_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_m_read_readvariableopZsavev2_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_m_read_readvariableopNsavev2_adam_shallow_sparse_model_5124_dense_48049_kernel_m_read_readvariableopLsavev2_adam_shallow_sparse_model_5124_dense_48049_bias_m_read_readvariableopNsavev2_adam_shallow_sparse_model_5124_dense_48048_kernel_v_read_readvariableopLsavev2_adam_shallow_sparse_model_5124_dense_48048_bias_v_read_readvariableop[savev2_adam_shallow_sparse_model_5124_batch_normalization_32124_gamma_v_read_readvariableopZsavev2_adam_shallow_sparse_model_5124_batch_normalization_32124_beta_v_read_readvariableopNsavev2_adam_shallow_sparse_model_5124_dense_48049_kernel_v_read_readvariableopLsavev2_adam_shallow_sparse_model_5124_dense_48049_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *,
dtypes"
 2	�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: Q

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: [
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 "!

identity_1Identity_1:output:0*�
_input_shapes�
�: :	�Z:Z:Z:Z:Z:Z:Z:: : : : : : : : : :	�Z:Z:Z:Z:Z::	�Z:Z:Z:Z:Z:: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:%!

_output_shapes
:	�Z: 

_output_shapes
:Z: 

_output_shapes
:Z: 

_output_shapes
:Z: 

_output_shapes
:Z: 

_output_shapes
:Z:$ 

_output_shapes

:Z: 

_output_shapes
::	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :%!

_output_shapes
:	�Z: 

_output_shapes
:Z: 

_output_shapes
:Z: 

_output_shapes
:Z:$ 

_output_shapes

:Z: 

_output_shapes
::%!

_output_shapes
:	�Z: 

_output_shapes
:Z: 

_output_shapes
:Z: 

_output_shapes
:Z:$ 

_output_shapes

:Z: 

_output_shapes
::

_output_shapes
: 
�

�
I__inference_dense_48049_layer_call_and_return_conditional_losses_57494956

inputs0
matmul_readvariableop_resource:Z-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:Z*
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
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������Z: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�
�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495134
input_1'
dense_48048_57495114:	�Z"
dense_48048_57495116:Z0
"batch_normalization_32124_57495119:Z0
"batch_normalization_32124_57495121:Z0
"batch_normalization_32124_57495123:Z0
"batch_normalization_32124_57495125:Z&
dense_48049_57495128:Z"
dense_48049_57495130:
identity��1batch_normalization_32124/StatefulPartitionedCall�#dense_48048/StatefulPartitionedCall�#dense_48049/StatefulPartitionedCall�
#dense_48048/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_48048_57495114dense_48048_57495116*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_48048_layer_call_and_return_conditional_losses_57494930�
1batch_normalization_32124/StatefulPartitionedCallStatefulPartitionedCall,dense_48048/StatefulPartitionedCall:output:0"batch_normalization_32124_57495119"batch_normalization_32124_57495121"batch_normalization_32124_57495123"batch_normalization_32124_57495125*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494902�
#dense_48049/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_32124/StatefulPartitionedCall:output:0dense_48049_57495128dense_48049_57495130*
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
GPU 2J 8� *R
fMRK
I__inference_dense_48049_layer_call_and_return_conditional_losses_57494956{
IdentityIdentity,dense_48049/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_32124/StatefulPartitionedCall$^dense_48048/StatefulPartitionedCall$^dense_48049/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_32124/StatefulPartitionedCall1batch_normalization_32124/StatefulPartitionedCall2J
#dense_48048/StatefulPartitionedCall#dense_48048/StatefulPartitionedCall2J
#dense_48049/StatefulPartitionedCall#dense_48049/StatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
<__inference_batch_normalization_32124_layer_call_fn_57495330

inputs
unknown:Z
	unknown_0:Z
	unknown_1:Z
	unknown_2:Z
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������Z*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494902o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������Z`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������Z: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs
�	
�
<__inference_shallow_sparse_model_5124_layer_call_fn_57495184
x
unknown:	�Z
	unknown_0:Z
	unknown_1:Z
	unknown_2:Z
	unknown_3:Z
	unknown_4:Z
	unknown_5:Z
	unknown_6:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57494963o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�,
�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495238
x=
*dense_48048_matmul_readvariableop_resource:	�Z9
+dense_48048_biasadd_readvariableop_resource:ZD
6batch_normalization_32124_cast_readvariableop_resource:ZF
8batch_normalization_32124_cast_1_readvariableop_resource:ZF
8batch_normalization_32124_cast_2_readvariableop_resource:ZF
8batch_normalization_32124_cast_3_readvariableop_resource:Z<
*dense_48049_matmul_readvariableop_resource:Z9
+dense_48049_biasadd_readvariableop_resource:
identity��-batch_normalization_32124/Cast/ReadVariableOp�/batch_normalization_32124/Cast_1/ReadVariableOp�/batch_normalization_32124/Cast_2/ReadVariableOp�/batch_normalization_32124/Cast_3/ReadVariableOp�"dense_48048/BiasAdd/ReadVariableOp�!dense_48048/MatMul/ReadVariableOp�"dense_48049/BiasAdd/ReadVariableOp�!dense_48049/MatMul/ReadVariableOp�
!dense_48048/MatMul/ReadVariableOpReadVariableOp*dense_48048_matmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0|
dense_48048/MatMulMatMulx)dense_48048/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
"dense_48048/BiasAdd/ReadVariableOpReadVariableOp+dense_48048_biasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0�
dense_48048/BiasAddBiasAdddense_48048/MatMul:product:0*dense_48048/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
-batch_normalization_32124/Cast/ReadVariableOpReadVariableOp6batch_normalization_32124_cast_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_32124/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_32124_cast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_32124/Cast_2/ReadVariableOpReadVariableOp8batch_normalization_32124_cast_2_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_32124/Cast_3/ReadVariableOpReadVariableOp8batch_normalization_32124_cast_3_readvariableop_resource*
_output_shapes
:Z*
dtype0n
)batch_normalization_32124/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_32124/batchnorm/addAddV27batch_normalization_32124/Cast_1/ReadVariableOp:value:02batch_normalization_32124/batchnorm/add/y:output:0*
T0*
_output_shapes
:Z�
)batch_normalization_32124/batchnorm/RsqrtRsqrt+batch_normalization_32124/batchnorm/add:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_32124/batchnorm/mulMul-batch_normalization_32124/batchnorm/Rsqrt:y:07batch_normalization_32124/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:Z�
)batch_normalization_32124/batchnorm/mul_1Muldense_48048/BiasAdd:output:0+batch_normalization_32124/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������Z�
)batch_normalization_32124/batchnorm/mul_2Mul5batch_normalization_32124/Cast/ReadVariableOp:value:0+batch_normalization_32124/batchnorm/mul:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_32124/batchnorm/subSub7batch_normalization_32124/Cast_2/ReadVariableOp:value:0-batch_normalization_32124/batchnorm/mul_2:z:0*
T0*
_output_shapes
:Z�
)batch_normalization_32124/batchnorm/add_1AddV2-batch_normalization_32124/batchnorm/mul_1:z:0+batch_normalization_32124/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Z�
!dense_48049/MatMul/ReadVariableOpReadVariableOp*dense_48049_matmul_readvariableop_resource*
_output_shapes

:Z*
dtype0�
dense_48049/MatMulMatMul-batch_normalization_32124/batchnorm/add_1:z:0)dense_48049/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
"dense_48049/BiasAdd/ReadVariableOpReadVariableOp+dense_48049_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_48049/BiasAddBiasAdddense_48049/MatMul:product:0*dense_48049/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������n
dense_48049/SigmoidSigmoiddense_48049/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
IdentityIdentitydense_48049/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp.^batch_normalization_32124/Cast/ReadVariableOp0^batch_normalization_32124/Cast_1/ReadVariableOp0^batch_normalization_32124/Cast_2/ReadVariableOp0^batch_normalization_32124/Cast_3/ReadVariableOp#^dense_48048/BiasAdd/ReadVariableOp"^dense_48048/MatMul/ReadVariableOp#^dense_48049/BiasAdd/ReadVariableOp"^dense_48049/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2^
-batch_normalization_32124/Cast/ReadVariableOp-batch_normalization_32124/Cast/ReadVariableOp2b
/batch_normalization_32124/Cast_1/ReadVariableOp/batch_normalization_32124/Cast_1/ReadVariableOp2b
/batch_normalization_32124/Cast_2/ReadVariableOp/batch_normalization_32124/Cast_2/ReadVariableOp2b
/batch_normalization_32124/Cast_3/ReadVariableOp/batch_normalization_32124/Cast_3/ReadVariableOp2H
"dense_48048/BiasAdd/ReadVariableOp"dense_48048/BiasAdd/ReadVariableOp2F
!dense_48048/MatMul/ReadVariableOp!dense_48048/MatMul/ReadVariableOp2H
"dense_48049/BiasAdd/ReadVariableOp"dense_48049/BiasAdd/ReadVariableOp2F
!dense_48049/MatMul/ReadVariableOp!dense_48049/MatMul/ReadVariableOp:K G
(
_output_shapes
:����������

_user_specified_namex
�
�
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57494855

inputs*
cast_readvariableop_resource:Z,
cast_1_readvariableop_resource:Z,
cast_2_readvariableop_resource:Z,
cast_3_readvariableop_resource:Z
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:Z*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:Z*
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
:ZP
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:Zm
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:Zc
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������Zk
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:Zm
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:Zr
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Zb
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������Z�
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������Z: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:O K
'
_output_shapes
:���������Z
 
_user_specified_nameinputs"�	L
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
StatefulPartitionedCall:0���������tensorflow/serving/predict:ā
�
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	embed
		norm1

out
	optimizer

signatures"
_tf_keras_model
X
0
1
2
3
4
5
6
7"
trackable_list_wrapper
J
0
1
2
3
4
5"
trackable_list_wrapper
 "
trackable_list_wrapper
�
non_trainable_variables

layers
metrics
layer_regularization_losses
layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�
trace_0
trace_1
trace_2
trace_32�
<__inference_shallow_sparse_model_5124_layer_call_fn_57494982
<__inference_shallow_sparse_model_5124_layer_call_fn_57495184
<__inference_shallow_sparse_model_5124_layer_call_fn_57495205
<__inference_shallow_sparse_model_5124_layer_call_fn_57495088�
���
FullArgSpec
args�
jself
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
annotations� *
 ztrace_0ztrace_1ztrace_2ztrace_3
�
trace_0
trace_1
 trace_2
!trace_32�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495238
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495285
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495111
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495134�
���
FullArgSpec
args�
jself
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
annotations� *
 ztrace_0ztrace_1z trace_2z!trace_3
�B�
#__inference__wrapped_model_57494831input_1"�
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
"	variables
#trainable_variables
$regularization_losses
%	keras_api
&__call__
*'&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
(	variables
)trainable_variables
*regularization_losses
+	keras_api
,__call__
*-&call_and_return_all_conditional_losses
.axis
	gamma
beta
moving_mean
moving_variance"
_tf_keras_layer
�
/	variables
0trainable_variables
1regularization_losses
2	keras_api
3__call__
*4&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
5iter

6beta_1

7beta_2
	8decay
9learning_ratem]m^m_m`mambvcvdvevfvgvh"
	optimizer
,
:serving_default"
signature_map
?:=	�Z2,shallow_sparse_model_5124/dense_48048/kernel
8:6Z2*shallow_sparse_model_5124/dense_48048/bias
G:EZ29shallow_sparse_model_5124/batch_normalization_32124/gamma
F:DZ28shallow_sparse_model_5124/batch_normalization_32124/beta
O:MZ (2?shallow_sparse_model_5124/batch_normalization_32124/moving_mean
S:QZ (2Cshallow_sparse_model_5124/batch_normalization_32124/moving_variance
>:<Z2,shallow_sparse_model_5124/dense_48049/kernel
8:62*shallow_sparse_model_5124/dense_48049/bias
.
0
1"
trackable_list_wrapper
5
0
	1

2"
trackable_list_wrapper
.
;0
<1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
<__inference_shallow_sparse_model_5124_layer_call_fn_57494982input_1"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
<__inference_shallow_sparse_model_5124_layer_call_fn_57495184x"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
<__inference_shallow_sparse_model_5124_layer_call_fn_57495205x"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
<__inference_shallow_sparse_model_5124_layer_call_fn_57495088input_1"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495238x"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495285x"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495111input_1"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495134input_1"�
���
FullArgSpec
args�
jself
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
annotations� *
 
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
=non_trainable_variables

>layers
?metrics
@layer_regularization_losses
Alayer_metrics
"	variables
#trainable_variables
$regularization_losses
&__call__
*'&call_and_return_all_conditional_losses
&'"call_and_return_conditional_losses"
_generic_user_object
�
Btrace_02�
.__inference_dense_48048_layer_call_fn_57495294�
���
FullArgSpec
args�
jself
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
annotations� *
 zBtrace_0
�
Ctrace_02�
I__inference_dense_48048_layer_call_and_return_conditional_losses_57495304�
���
FullArgSpec
args�
jself
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
annotations� *
 zCtrace_0
<
0
1
2
3"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
Dnon_trainable_variables

Elayers
Fmetrics
Glayer_regularization_losses
Hlayer_metrics
(	variables
)trainable_variables
*regularization_losses
,__call__
*-&call_and_return_all_conditional_losses
&-"call_and_return_conditional_losses"
_generic_user_object
�
Itrace_0
Jtrace_12�
<__inference_batch_normalization_32124_layer_call_fn_57495317
<__inference_batch_normalization_32124_layer_call_fn_57495330�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zItrace_0zJtrace_1
�
Ktrace_0
Ltrace_12�
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495350
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495384�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zKtrace_0zLtrace_1
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
Mnon_trainable_variables

Nlayers
Ometrics
Player_regularization_losses
Qlayer_metrics
/	variables
0trainable_variables
1regularization_losses
3__call__
*4&call_and_return_all_conditional_losses
&4"call_and_return_conditional_losses"
_generic_user_object
�
Rtrace_02�
.__inference_dense_48049_layer_call_fn_57495393�
���
FullArgSpec
args�
jself
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
annotations� *
 zRtrace_0
�
Strace_02�
I__inference_dense_48049_layer_call_and_return_conditional_losses_57495404�
���
FullArgSpec
args�
jself
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
annotations� *
 zStrace_0
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
�B�
&__inference_signature_wrapper_57495163input_1"�
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
N
T	variables
U	keras_api
	Vtotal
	Wcount"
_tf_keras_metric
^
X	variables
Y	keras_api
	Ztotal
	[count
\
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
.__inference_dense_48048_layer_call_fn_57495294inputs"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
I__inference_dense_48048_layer_call_and_return_conditional_losses_57495304inputs"�
���
FullArgSpec
args�
jself
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
annotations� *
 
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
<__inference_batch_normalization_32124_layer_call_fn_57495317inputs"�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
<__inference_batch_normalization_32124_layer_call_fn_57495330inputs"�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495350inputs"�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495384inputs"�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
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
.__inference_dense_48049_layer_call_fn_57495393inputs"�
���
FullArgSpec
args�
jself
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
annotations� *
 
�B�
I__inference_dense_48049_layer_call_and_return_conditional_losses_57495404inputs"�
���
FullArgSpec
args�
jself
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
annotations� *
 
.
V0
W1"
trackable_list_wrapper
-
T	variables"
_generic_user_object
:  (2total
:  (2count
.
Z0
[1"
trackable_list_wrapper
-
X	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
D:B	�Z23Adam/shallow_sparse_model_5124/dense_48048/kernel/m
=:;Z21Adam/shallow_sparse_model_5124/dense_48048/bias/m
L:JZ2@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/m
K:IZ2?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/m
C:AZ23Adam/shallow_sparse_model_5124/dense_48049/kernel/m
=:;21Adam/shallow_sparse_model_5124/dense_48049/bias/m
D:B	�Z23Adam/shallow_sparse_model_5124/dense_48048/kernel/v
=:;Z21Adam/shallow_sparse_model_5124/dense_48048/bias/v
L:JZ2@Adam/shallow_sparse_model_5124/batch_normalization_32124/gamma/v
K:IZ2?Adam/shallow_sparse_model_5124/batch_normalization_32124/beta/v
C:AZ23Adam/shallow_sparse_model_5124/dense_48049/kernel/v
=:;21Adam/shallow_sparse_model_5124/dense_48049/bias/v�
#__inference__wrapped_model_57494831r1�.
'�$
"�
input_1����������
� "3�0
.
output_1"�
output_1����������
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495350b3�0
)�&
 �
inputs���������Z
p 
� "%�"
�
0���������Z
� �
W__inference_batch_normalization_32124_layer_call_and_return_conditional_losses_57495384b3�0
)�&
 �
inputs���������Z
p
� "%�"
�
0���������Z
� �
<__inference_batch_normalization_32124_layer_call_fn_57495317U3�0
)�&
 �
inputs���������Z
p 
� "����������Z�
<__inference_batch_normalization_32124_layer_call_fn_57495330U3�0
)�&
 �
inputs���������Z
p
� "����������Z�
I__inference_dense_48048_layer_call_and_return_conditional_losses_57495304]0�-
&�#
!�
inputs����������
� "%�"
�
0���������Z
� �
.__inference_dense_48048_layer_call_fn_57495294P0�-
&�#
!�
inputs����������
� "����������Z�
I__inference_dense_48049_layer_call_and_return_conditional_losses_57495404\/�,
%�"
 �
inputs���������Z
� "%�"
�
0���������
� �
.__inference_dense_48049_layer_call_fn_57495393O/�,
%�"
 �
inputs���������Z
� "�����������
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495111tA�>
'�$
"�
input_1����������
�

trainingp "%�"
�
0���������
� �
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495134tA�>
'�$
"�
input_1����������
�

trainingp"%�"
�
0���������
� �
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495238n;�8
!�
�
x����������
�

trainingp "%�"
�
0���������
� �
W__inference_shallow_sparse_model_5124_layer_call_and_return_conditional_losses_57495285n;�8
!�
�
x����������
�

trainingp"%�"
�
0���������
� �
<__inference_shallow_sparse_model_5124_layer_call_fn_57494982gA�>
'�$
"�
input_1����������
�

trainingp "�����������
<__inference_shallow_sparse_model_5124_layer_call_fn_57495088gA�>
'�$
"�
input_1����������
�

trainingp"�����������
<__inference_shallow_sparse_model_5124_layer_call_fn_57495184a;�8
!�
�
x����������
�

trainingp "�����������
<__inference_shallow_sparse_model_5124_layer_call_fn_57495205a;�8
!�
�
x����������
�

trainingp"�����������
&__inference_signature_wrapper_57495163}<�9
� 
2�/
-
input_1"�
input_1����������"3�0
.
output_1"�
output_1���������