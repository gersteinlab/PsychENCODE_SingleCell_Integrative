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
2Adam/shallow_sparse_model_10528/dense_96657/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/shallow_sparse_model_10528/dense_96657/bias/v
�
FAdam/shallow_sparse_model_10528/dense_96657/bias/v/Read/ReadVariableOpReadVariableOp2Adam/shallow_sparse_model_10528/dense_96657/bias/v*
_output_shapes
:*
dtype0
�
4Adam/shallow_sparse_model_10528/dense_96657/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:Z*E
shared_name64Adam/shallow_sparse_model_10528/dense_96657/kernel/v
�
HAdam/shallow_sparse_model_10528/dense_96657/kernel/v/Read/ReadVariableOpReadVariableOp4Adam/shallow_sparse_model_10528/dense_96657/kernel/v*
_output_shapes

:Z*
dtype0
�
@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*Q
shared_nameB@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/v
�
TAdam/shallow_sparse_model_10528/batch_normalization_64528/beta/v/Read/ReadVariableOpReadVariableOp@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/v*
_output_shapes
:Z*
dtype0
�
AAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*R
shared_nameCAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/v
�
UAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/v/Read/ReadVariableOpReadVariableOpAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/v*
_output_shapes
:Z*
dtype0
�
2Adam/shallow_sparse_model_10528/dense_96656/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*C
shared_name42Adam/shallow_sparse_model_10528/dense_96656/bias/v
�
FAdam/shallow_sparse_model_10528/dense_96656/bias/v/Read/ReadVariableOpReadVariableOp2Adam/shallow_sparse_model_10528/dense_96656/bias/v*
_output_shapes
:Z*
dtype0
�
4Adam/shallow_sparse_model_10528/dense_96656/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�Z*E
shared_name64Adam/shallow_sparse_model_10528/dense_96656/kernel/v
�
HAdam/shallow_sparse_model_10528/dense_96656/kernel/v/Read/ReadVariableOpReadVariableOp4Adam/shallow_sparse_model_10528/dense_96656/kernel/v*
_output_shapes
:	�Z*
dtype0
�
2Adam/shallow_sparse_model_10528/dense_96657/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/shallow_sparse_model_10528/dense_96657/bias/m
�
FAdam/shallow_sparse_model_10528/dense_96657/bias/m/Read/ReadVariableOpReadVariableOp2Adam/shallow_sparse_model_10528/dense_96657/bias/m*
_output_shapes
:*
dtype0
�
4Adam/shallow_sparse_model_10528/dense_96657/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:Z*E
shared_name64Adam/shallow_sparse_model_10528/dense_96657/kernel/m
�
HAdam/shallow_sparse_model_10528/dense_96657/kernel/m/Read/ReadVariableOpReadVariableOp4Adam/shallow_sparse_model_10528/dense_96657/kernel/m*
_output_shapes

:Z*
dtype0
�
@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*Q
shared_nameB@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/m
�
TAdam/shallow_sparse_model_10528/batch_normalization_64528/beta/m/Read/ReadVariableOpReadVariableOp@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/m*
_output_shapes
:Z*
dtype0
�
AAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*R
shared_nameCAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/m
�
UAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/m/Read/ReadVariableOpReadVariableOpAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/m*
_output_shapes
:Z*
dtype0
�
2Adam/shallow_sparse_model_10528/dense_96656/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*C
shared_name42Adam/shallow_sparse_model_10528/dense_96656/bias/m
�
FAdam/shallow_sparse_model_10528/dense_96656/bias/m/Read/ReadVariableOpReadVariableOp2Adam/shallow_sparse_model_10528/dense_96656/bias/m*
_output_shapes
:Z*
dtype0
�
4Adam/shallow_sparse_model_10528/dense_96656/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�Z*E
shared_name64Adam/shallow_sparse_model_10528/dense_96656/kernel/m
�
HAdam/shallow_sparse_model_10528/dense_96656/kernel/m/Read/ReadVariableOpReadVariableOp4Adam/shallow_sparse_model_10528/dense_96656/kernel/m*
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
+shallow_sparse_model_10528/dense_96657/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*<
shared_name-+shallow_sparse_model_10528/dense_96657/bias
�
?shallow_sparse_model_10528/dense_96657/bias/Read/ReadVariableOpReadVariableOp+shallow_sparse_model_10528/dense_96657/bias*
_output_shapes
:*
dtype0
�
-shallow_sparse_model_10528/dense_96657/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:Z*>
shared_name/-shallow_sparse_model_10528/dense_96657/kernel
�
Ashallow_sparse_model_10528/dense_96657/kernel/Read/ReadVariableOpReadVariableOp-shallow_sparse_model_10528/dense_96657/kernel*
_output_shapes

:Z*
dtype0
�
Dshallow_sparse_model_10528/batch_normalization_64528/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*U
shared_nameFDshallow_sparse_model_10528/batch_normalization_64528/moving_variance
�
Xshallow_sparse_model_10528/batch_normalization_64528/moving_variance/Read/ReadVariableOpReadVariableOpDshallow_sparse_model_10528/batch_normalization_64528/moving_variance*
_output_shapes
:Z*
dtype0
�
@shallow_sparse_model_10528/batch_normalization_64528/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*Q
shared_nameB@shallow_sparse_model_10528/batch_normalization_64528/moving_mean
�
Tshallow_sparse_model_10528/batch_normalization_64528/moving_mean/Read/ReadVariableOpReadVariableOp@shallow_sparse_model_10528/batch_normalization_64528/moving_mean*
_output_shapes
:Z*
dtype0
�
9shallow_sparse_model_10528/batch_normalization_64528/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*J
shared_name;9shallow_sparse_model_10528/batch_normalization_64528/beta
�
Mshallow_sparse_model_10528/batch_normalization_64528/beta/Read/ReadVariableOpReadVariableOp9shallow_sparse_model_10528/batch_normalization_64528/beta*
_output_shapes
:Z*
dtype0
�
:shallow_sparse_model_10528/batch_normalization_64528/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*K
shared_name<:shallow_sparse_model_10528/batch_normalization_64528/gamma
�
Nshallow_sparse_model_10528/batch_normalization_64528/gamma/Read/ReadVariableOpReadVariableOp:shallow_sparse_model_10528/batch_normalization_64528/gamma*
_output_shapes
:Z*
dtype0
�
+shallow_sparse_model_10528/dense_96656/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:Z*<
shared_name-+shallow_sparse_model_10528/dense_96656/bias
�
?shallow_sparse_model_10528/dense_96656/bias/Read/ReadVariableOpReadVariableOp+shallow_sparse_model_10528/dense_96656/bias*
_output_shapes
:Z*
dtype0
�
-shallow_sparse_model_10528/dense_96656/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�Z*>
shared_name/-shallow_sparse_model_10528/dense_96656/kernel
�
Ashallow_sparse_model_10528/dense_96656/kernel/Read/ReadVariableOpReadVariableOp-shallow_sparse_model_10528/dense_96656/kernel*
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
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1-shallow_sparse_model_10528/dense_96656/kernel+shallow_sparse_model_10528/dense_96656/bias@shallow_sparse_model_10528/batch_normalization_64528/moving_meanDshallow_sparse_model_10528/batch_normalization_64528/moving_variance9shallow_sparse_model_10528/batch_normalization_64528/beta:shallow_sparse_model_10528/batch_normalization_64528/gamma-shallow_sparse_model_10528/dense_96657/kernel+shallow_sparse_model_10528/dense_96657/bias*
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
GPU 2J 8� *0
f+R)
'__inference_signature_wrapper_115899219

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
mg
VARIABLE_VALUE-shallow_sparse_model_10528/dense_96656/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
ke
VARIABLE_VALUE+shallow_sparse_model_10528/dense_96656/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
zt
VARIABLE_VALUE:shallow_sparse_model_10528/batch_normalization_64528/gamma&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
ys
VARIABLE_VALUE9shallow_sparse_model_10528/batch_normalization_64528/beta&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUE@shallow_sparse_model_10528/batch_normalization_64528/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUEDshallow_sparse_model_10528/batch_normalization_64528/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
mg
VARIABLE_VALUE-shallow_sparse_model_10528/dense_96657/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE*
ke
VARIABLE_VALUE+shallow_sparse_model_10528/dense_96657/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE*
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
VARIABLE_VALUE4Adam/shallow_sparse_model_10528/dense_96656/kernel/mBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE2Adam/shallow_sparse_model_10528/dense_96656/bias/mBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUEAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/mBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/mBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE4Adam/shallow_sparse_model_10528/dense_96657/kernel/mBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE2Adam/shallow_sparse_model_10528/dense_96657/bias/mBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE4Adam/shallow_sparse_model_10528/dense_96656/kernel/vBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE2Adam/shallow_sparse_model_10528/dense_96656/bias/vBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUEAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/vBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/vBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE4Adam/shallow_sparse_model_10528/dense_96657/kernel/vBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE2Adam/shallow_sparse_model_10528/dense_96657/bias/vBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenameAshallow_sparse_model_10528/dense_96656/kernel/Read/ReadVariableOp?shallow_sparse_model_10528/dense_96656/bias/Read/ReadVariableOpNshallow_sparse_model_10528/batch_normalization_64528/gamma/Read/ReadVariableOpMshallow_sparse_model_10528/batch_normalization_64528/beta/Read/ReadVariableOpTshallow_sparse_model_10528/batch_normalization_64528/moving_mean/Read/ReadVariableOpXshallow_sparse_model_10528/batch_normalization_64528/moving_variance/Read/ReadVariableOpAshallow_sparse_model_10528/dense_96657/kernel/Read/ReadVariableOp?shallow_sparse_model_10528/dense_96657/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpHAdam/shallow_sparse_model_10528/dense_96656/kernel/m/Read/ReadVariableOpFAdam/shallow_sparse_model_10528/dense_96656/bias/m/Read/ReadVariableOpUAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/m/Read/ReadVariableOpTAdam/shallow_sparse_model_10528/batch_normalization_64528/beta/m/Read/ReadVariableOpHAdam/shallow_sparse_model_10528/dense_96657/kernel/m/Read/ReadVariableOpFAdam/shallow_sparse_model_10528/dense_96657/bias/m/Read/ReadVariableOpHAdam/shallow_sparse_model_10528/dense_96656/kernel/v/Read/ReadVariableOpFAdam/shallow_sparse_model_10528/dense_96656/bias/v/Read/ReadVariableOpUAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/v/Read/ReadVariableOpTAdam/shallow_sparse_model_10528/batch_normalization_64528/beta/v/Read/ReadVariableOpHAdam/shallow_sparse_model_10528/dense_96657/kernel/v/Read/ReadVariableOpFAdam/shallow_sparse_model_10528/dense_96657/bias/v/Read/ReadVariableOpConst**
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
GPU 2J 8� *+
f&R$
"__inference__traced_save_115899570
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename-shallow_sparse_model_10528/dense_96656/kernel+shallow_sparse_model_10528/dense_96656/bias:shallow_sparse_model_10528/batch_normalization_64528/gamma9shallow_sparse_model_10528/batch_normalization_64528/beta@shallow_sparse_model_10528/batch_normalization_64528/moving_meanDshallow_sparse_model_10528/batch_normalization_64528/moving_variance-shallow_sparse_model_10528/dense_96657/kernel+shallow_sparse_model_10528/dense_96657/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_ratetotal_1count_1totalcount4Adam/shallow_sparse_model_10528/dense_96656/kernel/m2Adam/shallow_sparse_model_10528/dense_96656/bias/mAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/m@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/m4Adam/shallow_sparse_model_10528/dense_96657/kernel/m2Adam/shallow_sparse_model_10528/dense_96657/bias/m4Adam/shallow_sparse_model_10528/dense_96656/kernel/v2Adam/shallow_sparse_model_10528/dense_96656/bias/vAAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/v@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/v4Adam/shallow_sparse_model_10528/dense_96657/kernel/v2Adam/shallow_sparse_model_10528/dense_96657/bias/v*)
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
GPU 2J 8� *.
f)R'
%__inference__traced_restore_115899667��
�	
�
'__inference_signature_wrapper_115899219
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
GPU 2J 8� *-
f(R&
$__inference__wrapped_model_115898887o
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
�G
�
"__inference__traced_save_115899570
file_prefixL
Hsavev2_shallow_sparse_model_10528_dense_96656_kernel_read_readvariableopJ
Fsavev2_shallow_sparse_model_10528_dense_96656_bias_read_readvariableopY
Usavev2_shallow_sparse_model_10528_batch_normalization_64528_gamma_read_readvariableopX
Tsavev2_shallow_sparse_model_10528_batch_normalization_64528_beta_read_readvariableop_
[savev2_shallow_sparse_model_10528_batch_normalization_64528_moving_mean_read_readvariableopc
_savev2_shallow_sparse_model_10528_batch_normalization_64528_moving_variance_read_readvariableopL
Hsavev2_shallow_sparse_model_10528_dense_96657_kernel_read_readvariableopJ
Fsavev2_shallow_sparse_model_10528_dense_96657_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableopS
Osavev2_adam_shallow_sparse_model_10528_dense_96656_kernel_m_read_readvariableopQ
Msavev2_adam_shallow_sparse_model_10528_dense_96656_bias_m_read_readvariableop`
\savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_m_read_readvariableop_
[savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_m_read_readvariableopS
Osavev2_adam_shallow_sparse_model_10528_dense_96657_kernel_m_read_readvariableopQ
Msavev2_adam_shallow_sparse_model_10528_dense_96657_bias_m_read_readvariableopS
Osavev2_adam_shallow_sparse_model_10528_dense_96656_kernel_v_read_readvariableopQ
Msavev2_adam_shallow_sparse_model_10528_dense_96656_bias_v_read_readvariableop`
\savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_v_read_readvariableop_
[savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_v_read_readvariableopS
Osavev2_adam_shallow_sparse_model_10528_dense_96657_kernel_v_read_readvariableopQ
Msavev2_adam_shallow_sparse_model_10528_dense_96657_bias_v_read_readvariableop
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
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0Hsavev2_shallow_sparse_model_10528_dense_96656_kernel_read_readvariableopFsavev2_shallow_sparse_model_10528_dense_96656_bias_read_readvariableopUsavev2_shallow_sparse_model_10528_batch_normalization_64528_gamma_read_readvariableopTsavev2_shallow_sparse_model_10528_batch_normalization_64528_beta_read_readvariableop[savev2_shallow_sparse_model_10528_batch_normalization_64528_moving_mean_read_readvariableop_savev2_shallow_sparse_model_10528_batch_normalization_64528_moving_variance_read_readvariableopHsavev2_shallow_sparse_model_10528_dense_96657_kernel_read_readvariableopFsavev2_shallow_sparse_model_10528_dense_96657_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopOsavev2_adam_shallow_sparse_model_10528_dense_96656_kernel_m_read_readvariableopMsavev2_adam_shallow_sparse_model_10528_dense_96656_bias_m_read_readvariableop\savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_m_read_readvariableop[savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_m_read_readvariableopOsavev2_adam_shallow_sparse_model_10528_dense_96657_kernel_m_read_readvariableopMsavev2_adam_shallow_sparse_model_10528_dense_96657_bias_m_read_readvariableopOsavev2_adam_shallow_sparse_model_10528_dense_96656_kernel_v_read_readvariableopMsavev2_adam_shallow_sparse_model_10528_dense_96656_bias_v_read_readvariableop\savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_v_read_readvariableop[savev2_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_v_read_readvariableopOsavev2_adam_shallow_sparse_model_10528_dense_96657_kernel_v_read_readvariableopMsavev2_adam_shallow_sparse_model_10528_dense_96657_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
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
�
�
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898911

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
�?
�

$__inference__wrapped_model_115898887
input_1X
Eshallow_sparse_model_10528_dense_96656_matmul_readvariableop_resource:	�ZT
Fshallow_sparse_model_10528_dense_96656_biasadd_readvariableop_resource:Z_
Qshallow_sparse_model_10528_batch_normalization_64528_cast_readvariableop_resource:Za
Sshallow_sparse_model_10528_batch_normalization_64528_cast_1_readvariableop_resource:Za
Sshallow_sparse_model_10528_batch_normalization_64528_cast_2_readvariableop_resource:Za
Sshallow_sparse_model_10528_batch_normalization_64528_cast_3_readvariableop_resource:ZW
Eshallow_sparse_model_10528_dense_96657_matmul_readvariableop_resource:ZT
Fshallow_sparse_model_10528_dense_96657_biasadd_readvariableop_resource:
identity��Hshallow_sparse_model_10528/batch_normalization_64528/Cast/ReadVariableOp�Jshallow_sparse_model_10528/batch_normalization_64528/Cast_1/ReadVariableOp�Jshallow_sparse_model_10528/batch_normalization_64528/Cast_2/ReadVariableOp�Jshallow_sparse_model_10528/batch_normalization_64528/Cast_3/ReadVariableOp�=shallow_sparse_model_10528/dense_96656/BiasAdd/ReadVariableOp�<shallow_sparse_model_10528/dense_96656/MatMul/ReadVariableOp�=shallow_sparse_model_10528/dense_96657/BiasAdd/ReadVariableOp�<shallow_sparse_model_10528/dense_96657/MatMul/ReadVariableOp�
<shallow_sparse_model_10528/dense_96656/MatMul/ReadVariableOpReadVariableOpEshallow_sparse_model_10528_dense_96656_matmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0�
-shallow_sparse_model_10528/dense_96656/MatMulMatMulinput_1Dshallow_sparse_model_10528/dense_96656/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
=shallow_sparse_model_10528/dense_96656/BiasAdd/ReadVariableOpReadVariableOpFshallow_sparse_model_10528_dense_96656_biasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0�
.shallow_sparse_model_10528/dense_96656/BiasAddBiasAdd7shallow_sparse_model_10528/dense_96656/MatMul:product:0Eshallow_sparse_model_10528/dense_96656/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
Hshallow_sparse_model_10528/batch_normalization_64528/Cast/ReadVariableOpReadVariableOpQshallow_sparse_model_10528_batch_normalization_64528_cast_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Jshallow_sparse_model_10528/batch_normalization_64528/Cast_1/ReadVariableOpReadVariableOpSshallow_sparse_model_10528_batch_normalization_64528_cast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Jshallow_sparse_model_10528/batch_normalization_64528/Cast_2/ReadVariableOpReadVariableOpSshallow_sparse_model_10528_batch_normalization_64528_cast_2_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Jshallow_sparse_model_10528/batch_normalization_64528/Cast_3/ReadVariableOpReadVariableOpSshallow_sparse_model_10528_batch_normalization_64528_cast_3_readvariableop_resource*
_output_shapes
:Z*
dtype0�
Dshallow_sparse_model_10528/batch_normalization_64528/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
Bshallow_sparse_model_10528/batch_normalization_64528/batchnorm/addAddV2Rshallow_sparse_model_10528/batch_normalization_64528/Cast_1/ReadVariableOp:value:0Mshallow_sparse_model_10528/batch_normalization_64528/batchnorm/add/y:output:0*
T0*
_output_shapes
:Z�
Dshallow_sparse_model_10528/batch_normalization_64528/batchnorm/RsqrtRsqrtFshallow_sparse_model_10528/batch_normalization_64528/batchnorm/add:z:0*
T0*
_output_shapes
:Z�
Bshallow_sparse_model_10528/batch_normalization_64528/batchnorm/mulMulHshallow_sparse_model_10528/batch_normalization_64528/batchnorm/Rsqrt:y:0Rshallow_sparse_model_10528/batch_normalization_64528/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:Z�
Dshallow_sparse_model_10528/batch_normalization_64528/batchnorm/mul_1Mul7shallow_sparse_model_10528/dense_96656/BiasAdd:output:0Fshallow_sparse_model_10528/batch_normalization_64528/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������Z�
Dshallow_sparse_model_10528/batch_normalization_64528/batchnorm/mul_2MulPshallow_sparse_model_10528/batch_normalization_64528/Cast/ReadVariableOp:value:0Fshallow_sparse_model_10528/batch_normalization_64528/batchnorm/mul:z:0*
T0*
_output_shapes
:Z�
Bshallow_sparse_model_10528/batch_normalization_64528/batchnorm/subSubRshallow_sparse_model_10528/batch_normalization_64528/Cast_2/ReadVariableOp:value:0Hshallow_sparse_model_10528/batch_normalization_64528/batchnorm/mul_2:z:0*
T0*
_output_shapes
:Z�
Dshallow_sparse_model_10528/batch_normalization_64528/batchnorm/add_1AddV2Hshallow_sparse_model_10528/batch_normalization_64528/batchnorm/mul_1:z:0Fshallow_sparse_model_10528/batch_normalization_64528/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Z�
<shallow_sparse_model_10528/dense_96657/MatMul/ReadVariableOpReadVariableOpEshallow_sparse_model_10528_dense_96657_matmul_readvariableop_resource*
_output_shapes

:Z*
dtype0�
-shallow_sparse_model_10528/dense_96657/MatMulMatMulHshallow_sparse_model_10528/batch_normalization_64528/batchnorm/add_1:z:0Dshallow_sparse_model_10528/dense_96657/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
=shallow_sparse_model_10528/dense_96657/BiasAdd/ReadVariableOpReadVariableOpFshallow_sparse_model_10528_dense_96657_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
.shallow_sparse_model_10528/dense_96657/BiasAddBiasAdd7shallow_sparse_model_10528/dense_96657/MatMul:product:0Eshallow_sparse_model_10528/dense_96657/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.shallow_sparse_model_10528/dense_96657/SigmoidSigmoid7shallow_sparse_model_10528/dense_96657/BiasAdd:output:0*
T0*'
_output_shapes
:����������
IdentityIdentity2shallow_sparse_model_10528/dense_96657/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOpI^shallow_sparse_model_10528/batch_normalization_64528/Cast/ReadVariableOpK^shallow_sparse_model_10528/batch_normalization_64528/Cast_1/ReadVariableOpK^shallow_sparse_model_10528/batch_normalization_64528/Cast_2/ReadVariableOpK^shallow_sparse_model_10528/batch_normalization_64528/Cast_3/ReadVariableOp>^shallow_sparse_model_10528/dense_96656/BiasAdd/ReadVariableOp=^shallow_sparse_model_10528/dense_96656/MatMul/ReadVariableOp>^shallow_sparse_model_10528/dense_96657/BiasAdd/ReadVariableOp=^shallow_sparse_model_10528/dense_96657/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2�
Hshallow_sparse_model_10528/batch_normalization_64528/Cast/ReadVariableOpHshallow_sparse_model_10528/batch_normalization_64528/Cast/ReadVariableOp2�
Jshallow_sparse_model_10528/batch_normalization_64528/Cast_1/ReadVariableOpJshallow_sparse_model_10528/batch_normalization_64528/Cast_1/ReadVariableOp2�
Jshallow_sparse_model_10528/batch_normalization_64528/Cast_2/ReadVariableOpJshallow_sparse_model_10528/batch_normalization_64528/Cast_2/ReadVariableOp2�
Jshallow_sparse_model_10528/batch_normalization_64528/Cast_3/ReadVariableOpJshallow_sparse_model_10528/batch_normalization_64528/Cast_3/ReadVariableOp2~
=shallow_sparse_model_10528/dense_96656/BiasAdd/ReadVariableOp=shallow_sparse_model_10528/dense_96656/BiasAdd/ReadVariableOp2|
<shallow_sparse_model_10528/dense_96656/MatMul/ReadVariableOp<shallow_sparse_model_10528/dense_96656/MatMul/ReadVariableOp2~
=shallow_sparse_model_10528/dense_96657/BiasAdd/ReadVariableOp=shallow_sparse_model_10528/dense_96657/BiasAdd/ReadVariableOp2|
<shallow_sparse_model_10528/dense_96657/MatMul/ReadVariableOp<shallow_sparse_model_10528/dense_96657/MatMul/ReadVariableOp:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�

�
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899460

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
�	
�
>__inference_shallow_sparse_model_10528_layer_call_fn_115899144
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
GPU 2J 8� *b
f]R[
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899104o
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
�
J__inference_dense_96656_layer_call_and_return_conditional_losses_115899360

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
�	
�
>__inference_shallow_sparse_model_10528_layer_call_fn_115899240
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
GPU 2J 8� *b
f]R[
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899019o
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
�
J__inference_dense_96656_layer_call_and_return_conditional_losses_115898986

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
�$
�
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898958

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
�
�
=__inference_batch_normalization_64528_layer_call_fn_115899386

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
GPU 2J 8� *a
f\RZ
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898958o
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

�
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899012

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
%__inference__traced_restore_115899667
file_prefixQ
>assignvariableop_shallow_sparse_model_10528_dense_96656_kernel:	�ZL
>assignvariableop_1_shallow_sparse_model_10528_dense_96656_bias:Z[
Massignvariableop_2_shallow_sparse_model_10528_batch_normalization_64528_gamma:ZZ
Lassignvariableop_3_shallow_sparse_model_10528_batch_normalization_64528_beta:Za
Sassignvariableop_4_shallow_sparse_model_10528_batch_normalization_64528_moving_mean:Ze
Wassignvariableop_5_shallow_sparse_model_10528_batch_normalization_64528_moving_variance:ZR
@assignvariableop_6_shallow_sparse_model_10528_dense_96657_kernel:ZL
>assignvariableop_7_shallow_sparse_model_10528_dense_96657_bias:&
assignvariableop_8_adam_iter:	 (
assignvariableop_9_adam_beta_1: )
assignvariableop_10_adam_beta_2: (
assignvariableop_11_adam_decay: 0
&assignvariableop_12_adam_learning_rate: %
assignvariableop_13_total_1: %
assignvariableop_14_count_1: #
assignvariableop_15_total: #
assignvariableop_16_count: [
Hassignvariableop_17_adam_shallow_sparse_model_10528_dense_96656_kernel_m:	�ZT
Fassignvariableop_18_adam_shallow_sparse_model_10528_dense_96656_bias_m:Zc
Uassignvariableop_19_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_m:Zb
Tassignvariableop_20_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_m:ZZ
Hassignvariableop_21_adam_shallow_sparse_model_10528_dense_96657_kernel_m:ZT
Fassignvariableop_22_adam_shallow_sparse_model_10528_dense_96657_bias_m:[
Hassignvariableop_23_adam_shallow_sparse_model_10528_dense_96656_kernel_v:	�ZT
Fassignvariableop_24_adam_shallow_sparse_model_10528_dense_96656_bias_v:Zc
Uassignvariableop_25_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_v:Zb
Tassignvariableop_26_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_v:ZZ
Hassignvariableop_27_adam_shallow_sparse_model_10528_dense_96657_kernel_v:ZT
Fassignvariableop_28_adam_shallow_sparse_model_10528_dense_96657_bias_v:
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
AssignVariableOpAssignVariableOp>assignvariableop_shallow_sparse_model_10528_dense_96656_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp>assignvariableop_1_shallow_sparse_model_10528_dense_96656_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOpMassignvariableop_2_shallow_sparse_model_10528_batch_normalization_64528_gammaIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpLassignvariableop_3_shallow_sparse_model_10528_batch_normalization_64528_betaIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpSassignvariableop_4_shallow_sparse_model_10528_batch_normalization_64528_moving_meanIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpWassignvariableop_5_shallow_sparse_model_10528_batch_normalization_64528_moving_varianceIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp@assignvariableop_6_shallow_sparse_model_10528_dense_96657_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp>assignvariableop_7_shallow_sparse_model_10528_dense_96657_biasIdentity_7:output:0"/device:CPU:0*
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
AssignVariableOp_17AssignVariableOpHassignvariableop_17_adam_shallow_sparse_model_10528_dense_96656_kernel_mIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOpFassignvariableop_18_adam_shallow_sparse_model_10528_dense_96656_bias_mIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOpUassignvariableop_19_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_mIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOpTassignvariableop_20_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_mIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOpHassignvariableop_21_adam_shallow_sparse_model_10528_dense_96657_kernel_mIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOpFassignvariableop_22_adam_shallow_sparse_model_10528_dense_96657_bias_mIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOpHassignvariableop_23_adam_shallow_sparse_model_10528_dense_96656_kernel_vIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOpFassignvariableop_24_adam_shallow_sparse_model_10528_dense_96656_bias_vIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOpUassignvariableop_25_adam_shallow_sparse_model_10528_batch_normalization_64528_gamma_vIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpTassignvariableop_26_adam_shallow_sparse_model_10528_batch_normalization_64528_beta_vIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOpHassignvariableop_27_adam_shallow_sparse_model_10528_dense_96657_kernel_vIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpFassignvariableop_28_adam_shallow_sparse_model_10528_dense_96657_bias_vIdentity_28:output:0"/device:CPU:0*
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
�
�
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899019
x(
dense_96656_115898987:	�Z#
dense_96656_115898989:Z1
#batch_normalization_64528_115898992:Z1
#batch_normalization_64528_115898994:Z1
#batch_normalization_64528_115898996:Z1
#batch_normalization_64528_115898998:Z'
dense_96657_115899013:Z#
dense_96657_115899015:
identity��1batch_normalization_64528/StatefulPartitionedCall�#dense_96656/StatefulPartitionedCall�#dense_96657/StatefulPartitionedCall�
#dense_96656/StatefulPartitionedCallStatefulPartitionedCallxdense_96656_115898987dense_96656_115898989*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96656_layer_call_and_return_conditional_losses_115898986�
1batch_normalization_64528/StatefulPartitionedCallStatefulPartitionedCall,dense_96656/StatefulPartitionedCall:output:0#batch_normalization_64528_115898992#batch_normalization_64528_115898994#batch_normalization_64528_115898996#batch_normalization_64528_115898998*
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
GPU 2J 8� *a
f\RZ
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898911�
#dense_96657/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_64528/StatefulPartitionedCall:output:0dense_96657_115899013dense_96657_115899015*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899012{
IdentityIdentity,dense_96657/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_64528/StatefulPartitionedCall$^dense_96656/StatefulPartitionedCall$^dense_96657/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_64528/StatefulPartitionedCall1batch_normalization_64528/StatefulPartitionedCall2J
#dense_96656/StatefulPartitionedCall#dense_96656/StatefulPartitionedCall2J
#dense_96657/StatefulPartitionedCall#dense_96657/StatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�$
�
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899440

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
�,
�
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899294
x=
*dense_96656_matmul_readvariableop_resource:	�Z9
+dense_96656_biasadd_readvariableop_resource:ZD
6batch_normalization_64528_cast_readvariableop_resource:ZF
8batch_normalization_64528_cast_1_readvariableop_resource:ZF
8batch_normalization_64528_cast_2_readvariableop_resource:ZF
8batch_normalization_64528_cast_3_readvariableop_resource:Z<
*dense_96657_matmul_readvariableop_resource:Z9
+dense_96657_biasadd_readvariableop_resource:
identity��-batch_normalization_64528/Cast/ReadVariableOp�/batch_normalization_64528/Cast_1/ReadVariableOp�/batch_normalization_64528/Cast_2/ReadVariableOp�/batch_normalization_64528/Cast_3/ReadVariableOp�"dense_96656/BiasAdd/ReadVariableOp�!dense_96656/MatMul/ReadVariableOp�"dense_96657/BiasAdd/ReadVariableOp�!dense_96657/MatMul/ReadVariableOp�
!dense_96656/MatMul/ReadVariableOpReadVariableOp*dense_96656_matmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0|
dense_96656/MatMulMatMulx)dense_96656/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
"dense_96656/BiasAdd/ReadVariableOpReadVariableOp+dense_96656_biasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0�
dense_96656/BiasAddBiasAdddense_96656/MatMul:product:0*dense_96656/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
-batch_normalization_64528/Cast/ReadVariableOpReadVariableOp6batch_normalization_64528_cast_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_64528/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_64528_cast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_64528/Cast_2/ReadVariableOpReadVariableOp8batch_normalization_64528_cast_2_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_64528/Cast_3/ReadVariableOpReadVariableOp8batch_normalization_64528_cast_3_readvariableop_resource*
_output_shapes
:Z*
dtype0n
)batch_normalization_64528/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_64528/batchnorm/addAddV27batch_normalization_64528/Cast_1/ReadVariableOp:value:02batch_normalization_64528/batchnorm/add/y:output:0*
T0*
_output_shapes
:Z�
)batch_normalization_64528/batchnorm/RsqrtRsqrt+batch_normalization_64528/batchnorm/add:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_64528/batchnorm/mulMul-batch_normalization_64528/batchnorm/Rsqrt:y:07batch_normalization_64528/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:Z�
)batch_normalization_64528/batchnorm/mul_1Muldense_96656/BiasAdd:output:0+batch_normalization_64528/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������Z�
)batch_normalization_64528/batchnorm/mul_2Mul5batch_normalization_64528/Cast/ReadVariableOp:value:0+batch_normalization_64528/batchnorm/mul:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_64528/batchnorm/subSub7batch_normalization_64528/Cast_2/ReadVariableOp:value:0-batch_normalization_64528/batchnorm/mul_2:z:0*
T0*
_output_shapes
:Z�
)batch_normalization_64528/batchnorm/add_1AddV2-batch_normalization_64528/batchnorm/mul_1:z:0+batch_normalization_64528/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Z�
!dense_96657/MatMul/ReadVariableOpReadVariableOp*dense_96657_matmul_readvariableop_resource*
_output_shapes

:Z*
dtype0�
dense_96657/MatMulMatMul-batch_normalization_64528/batchnorm/add_1:z:0)dense_96657/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
"dense_96657/BiasAdd/ReadVariableOpReadVariableOp+dense_96657_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_96657/BiasAddBiasAdddense_96657/MatMul:product:0*dense_96657/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������n
dense_96657/SigmoidSigmoiddense_96657/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
IdentityIdentitydense_96657/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp.^batch_normalization_64528/Cast/ReadVariableOp0^batch_normalization_64528/Cast_1/ReadVariableOp0^batch_normalization_64528/Cast_2/ReadVariableOp0^batch_normalization_64528/Cast_3/ReadVariableOp#^dense_96656/BiasAdd/ReadVariableOp"^dense_96656/MatMul/ReadVariableOp#^dense_96657/BiasAdd/ReadVariableOp"^dense_96657/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2^
-batch_normalization_64528/Cast/ReadVariableOp-batch_normalization_64528/Cast/ReadVariableOp2b
/batch_normalization_64528/Cast_1/ReadVariableOp/batch_normalization_64528/Cast_1/ReadVariableOp2b
/batch_normalization_64528/Cast_2/ReadVariableOp/batch_normalization_64528/Cast_2/ReadVariableOp2b
/batch_normalization_64528/Cast_3/ReadVariableOp/batch_normalization_64528/Cast_3/ReadVariableOp2H
"dense_96656/BiasAdd/ReadVariableOp"dense_96656/BiasAdd/ReadVariableOp2F
!dense_96656/MatMul/ReadVariableOp!dense_96656/MatMul/ReadVariableOp2H
"dense_96657/BiasAdd/ReadVariableOp"dense_96657/BiasAdd/ReadVariableOp2F
!dense_96657/MatMul/ReadVariableOp!dense_96657/MatMul/ReadVariableOp:K G
(
_output_shapes
:����������

_user_specified_namex
�
�
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899190
input_1(
dense_96656_115899170:	�Z#
dense_96656_115899172:Z1
#batch_normalization_64528_115899175:Z1
#batch_normalization_64528_115899177:Z1
#batch_normalization_64528_115899179:Z1
#batch_normalization_64528_115899181:Z'
dense_96657_115899184:Z#
dense_96657_115899186:
identity��1batch_normalization_64528/StatefulPartitionedCall�#dense_96656/StatefulPartitionedCall�#dense_96657/StatefulPartitionedCall�
#dense_96656/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_96656_115899170dense_96656_115899172*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96656_layer_call_and_return_conditional_losses_115898986�
1batch_normalization_64528/StatefulPartitionedCallStatefulPartitionedCall,dense_96656/StatefulPartitionedCall:output:0#batch_normalization_64528_115899175#batch_normalization_64528_115899177#batch_normalization_64528_115899179#batch_normalization_64528_115899181*
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
GPU 2J 8� *a
f\RZ
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898958�
#dense_96657/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_64528/StatefulPartitionedCall:output:0dense_96657_115899184dense_96657_115899186*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899012{
IdentityIdentity,dense_96657/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_64528/StatefulPartitionedCall$^dense_96656/StatefulPartitionedCall$^dense_96657/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_64528/StatefulPartitionedCall1batch_normalization_64528/StatefulPartitionedCall2J
#dense_96656/StatefulPartitionedCall#dense_96656/StatefulPartitionedCall2J
#dense_96657/StatefulPartitionedCall#dense_96657/StatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899167
input_1(
dense_96656_115899147:	�Z#
dense_96656_115899149:Z1
#batch_normalization_64528_115899152:Z1
#batch_normalization_64528_115899154:Z1
#batch_normalization_64528_115899156:Z1
#batch_normalization_64528_115899158:Z'
dense_96657_115899161:Z#
dense_96657_115899163:
identity��1batch_normalization_64528/StatefulPartitionedCall�#dense_96656/StatefulPartitionedCall�#dense_96657/StatefulPartitionedCall�
#dense_96656/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_96656_115899147dense_96656_115899149*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96656_layer_call_and_return_conditional_losses_115898986�
1batch_normalization_64528/StatefulPartitionedCallStatefulPartitionedCall,dense_96656/StatefulPartitionedCall:output:0#batch_normalization_64528_115899152#batch_normalization_64528_115899154#batch_normalization_64528_115899156#batch_normalization_64528_115899158*
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
GPU 2J 8� *a
f\RZ
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898911�
#dense_96657/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_64528/StatefulPartitionedCall:output:0dense_96657_115899161dense_96657_115899163*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899012{
IdentityIdentity,dense_96657/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_64528/StatefulPartitionedCall$^dense_96656/StatefulPartitionedCall$^dense_96657/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_64528/StatefulPartitionedCall1batch_normalization_64528/StatefulPartitionedCall2J
#dense_96656/StatefulPartitionedCall#dense_96656/StatefulPartitionedCall2J
#dense_96657/StatefulPartitionedCall#dense_96657/StatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
/__inference_dense_96657_layer_call_fn_115899449

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
GPU 2J 8� *S
fNRL
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899012o
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
�
�
=__inference_batch_normalization_64528_layer_call_fn_115899373

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
GPU 2J 8� *a
f\RZ
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898911o
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
�H
�
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899341
x=
*dense_96656_matmul_readvariableop_resource:	�Z9
+dense_96656_biasadd_readvariableop_resource:ZO
Abatch_normalization_64528_assignmovingavg_readvariableop_resource:ZQ
Cbatch_normalization_64528_assignmovingavg_1_readvariableop_resource:ZD
6batch_normalization_64528_cast_readvariableop_resource:ZF
8batch_normalization_64528_cast_1_readvariableop_resource:Z<
*dense_96657_matmul_readvariableop_resource:Z9
+dense_96657_biasadd_readvariableop_resource:
identity��)batch_normalization_64528/AssignMovingAvg�8batch_normalization_64528/AssignMovingAvg/ReadVariableOp�+batch_normalization_64528/AssignMovingAvg_1�:batch_normalization_64528/AssignMovingAvg_1/ReadVariableOp�-batch_normalization_64528/Cast/ReadVariableOp�/batch_normalization_64528/Cast_1/ReadVariableOp�"dense_96656/BiasAdd/ReadVariableOp�!dense_96656/MatMul/ReadVariableOp�"dense_96657/BiasAdd/ReadVariableOp�!dense_96657/MatMul/ReadVariableOp�
!dense_96656/MatMul/ReadVariableOpReadVariableOp*dense_96656_matmul_readvariableop_resource*
_output_shapes
:	�Z*
dtype0|
dense_96656/MatMulMatMulx)dense_96656/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
"dense_96656/BiasAdd/ReadVariableOpReadVariableOp+dense_96656_biasadd_readvariableop_resource*
_output_shapes
:Z*
dtype0�
dense_96656/BiasAddBiasAdddense_96656/MatMul:product:0*dense_96656/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������Z�
8batch_normalization_64528/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
&batch_normalization_64528/moments/meanMeandense_96656/BiasAdd:output:0Abatch_normalization_64528/moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(�
.batch_normalization_64528/moments/StopGradientStopGradient/batch_normalization_64528/moments/mean:output:0*
T0*
_output_shapes

:Z�
3batch_normalization_64528/moments/SquaredDifferenceSquaredDifferencedense_96656/BiasAdd:output:07batch_normalization_64528/moments/StopGradient:output:0*
T0*'
_output_shapes
:���������Z�
<batch_normalization_64528/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
*batch_normalization_64528/moments/varianceMean7batch_normalization_64528/moments/SquaredDifference:z:0Ebatch_normalization_64528/moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:Z*
	keep_dims(�
)batch_normalization_64528/moments/SqueezeSqueeze/batch_normalization_64528/moments/mean:output:0*
T0*
_output_shapes
:Z*
squeeze_dims
 �
+batch_normalization_64528/moments/Squeeze_1Squeeze3batch_normalization_64528/moments/variance:output:0*
T0*
_output_shapes
:Z*
squeeze_dims
 t
/batch_normalization_64528/AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
8batch_normalization_64528/AssignMovingAvg/ReadVariableOpReadVariableOpAbatch_normalization_64528_assignmovingavg_readvariableop_resource*
_output_shapes
:Z*
dtype0�
-batch_normalization_64528/AssignMovingAvg/subSub@batch_normalization_64528/AssignMovingAvg/ReadVariableOp:value:02batch_normalization_64528/moments/Squeeze:output:0*
T0*
_output_shapes
:Z�
-batch_normalization_64528/AssignMovingAvg/mulMul1batch_normalization_64528/AssignMovingAvg/sub:z:08batch_normalization_64528/AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:Z�
)batch_normalization_64528/AssignMovingAvgAssignSubVariableOpAbatch_normalization_64528_assignmovingavg_readvariableop_resource1batch_normalization_64528/AssignMovingAvg/mul:z:09^batch_normalization_64528/AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0v
1batch_normalization_64528/AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
:batch_normalization_64528/AssignMovingAvg_1/ReadVariableOpReadVariableOpCbatch_normalization_64528_assignmovingavg_1_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_64528/AssignMovingAvg_1/subSubBbatch_normalization_64528/AssignMovingAvg_1/ReadVariableOp:value:04batch_normalization_64528/moments/Squeeze_1:output:0*
T0*
_output_shapes
:Z�
/batch_normalization_64528/AssignMovingAvg_1/mulMul3batch_normalization_64528/AssignMovingAvg_1/sub:z:0:batch_normalization_64528/AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:Z�
+batch_normalization_64528/AssignMovingAvg_1AssignSubVariableOpCbatch_normalization_64528_assignmovingavg_1_readvariableop_resource3batch_normalization_64528/AssignMovingAvg_1/mul:z:0;^batch_normalization_64528/AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0�
-batch_normalization_64528/Cast/ReadVariableOpReadVariableOp6batch_normalization_64528_cast_readvariableop_resource*
_output_shapes
:Z*
dtype0�
/batch_normalization_64528/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_64528_cast_1_readvariableop_resource*
_output_shapes
:Z*
dtype0n
)batch_normalization_64528/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_64528/batchnorm/addAddV24batch_normalization_64528/moments/Squeeze_1:output:02batch_normalization_64528/batchnorm/add/y:output:0*
T0*
_output_shapes
:Z�
)batch_normalization_64528/batchnorm/RsqrtRsqrt+batch_normalization_64528/batchnorm/add:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_64528/batchnorm/mulMul-batch_normalization_64528/batchnorm/Rsqrt:y:07batch_normalization_64528/Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:Z�
)batch_normalization_64528/batchnorm/mul_1Muldense_96656/BiasAdd:output:0+batch_normalization_64528/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������Z�
)batch_normalization_64528/batchnorm/mul_2Mul2batch_normalization_64528/moments/Squeeze:output:0+batch_normalization_64528/batchnorm/mul:z:0*
T0*
_output_shapes
:Z�
'batch_normalization_64528/batchnorm/subSub5batch_normalization_64528/Cast/ReadVariableOp:value:0-batch_normalization_64528/batchnorm/mul_2:z:0*
T0*
_output_shapes
:Z�
)batch_normalization_64528/batchnorm/add_1AddV2-batch_normalization_64528/batchnorm/mul_1:z:0+batch_normalization_64528/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������Z�
!dense_96657/MatMul/ReadVariableOpReadVariableOp*dense_96657_matmul_readvariableop_resource*
_output_shapes

:Z*
dtype0�
dense_96657/MatMulMatMul-batch_normalization_64528/batchnorm/add_1:z:0)dense_96657/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
"dense_96657/BiasAdd/ReadVariableOpReadVariableOp+dense_96657_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_96657/BiasAddBiasAdddense_96657/MatMul:product:0*dense_96657/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������n
dense_96657/SigmoidSigmoiddense_96657/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
IdentityIdentitydense_96657/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp*^batch_normalization_64528/AssignMovingAvg9^batch_normalization_64528/AssignMovingAvg/ReadVariableOp,^batch_normalization_64528/AssignMovingAvg_1;^batch_normalization_64528/AssignMovingAvg_1/ReadVariableOp.^batch_normalization_64528/Cast/ReadVariableOp0^batch_normalization_64528/Cast_1/ReadVariableOp#^dense_96656/BiasAdd/ReadVariableOp"^dense_96656/MatMul/ReadVariableOp#^dense_96657/BiasAdd/ReadVariableOp"^dense_96657/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2V
)batch_normalization_64528/AssignMovingAvg)batch_normalization_64528/AssignMovingAvg2t
8batch_normalization_64528/AssignMovingAvg/ReadVariableOp8batch_normalization_64528/AssignMovingAvg/ReadVariableOp2Z
+batch_normalization_64528/AssignMovingAvg_1+batch_normalization_64528/AssignMovingAvg_12x
:batch_normalization_64528/AssignMovingAvg_1/ReadVariableOp:batch_normalization_64528/AssignMovingAvg_1/ReadVariableOp2^
-batch_normalization_64528/Cast/ReadVariableOp-batch_normalization_64528/Cast/ReadVariableOp2b
/batch_normalization_64528/Cast_1/ReadVariableOp/batch_normalization_64528/Cast_1/ReadVariableOp2H
"dense_96656/BiasAdd/ReadVariableOp"dense_96656/BiasAdd/ReadVariableOp2F
!dense_96656/MatMul/ReadVariableOp!dense_96656/MatMul/ReadVariableOp2H
"dense_96657/BiasAdd/ReadVariableOp"dense_96657/BiasAdd/ReadVariableOp2F
!dense_96657/MatMul/ReadVariableOp!dense_96657/MatMul/ReadVariableOp:K G
(
_output_shapes
:����������

_user_specified_namex
�	
�
>__inference_shallow_sparse_model_10528_layer_call_fn_115899261
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
GPU 2J 8� *b
f]R[
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899104o
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
�
�
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899104
x(
dense_96656_115899084:	�Z#
dense_96656_115899086:Z1
#batch_normalization_64528_115899089:Z1
#batch_normalization_64528_115899091:Z1
#batch_normalization_64528_115899093:Z1
#batch_normalization_64528_115899095:Z'
dense_96657_115899098:Z#
dense_96657_115899100:
identity��1batch_normalization_64528/StatefulPartitionedCall�#dense_96656/StatefulPartitionedCall�#dense_96657/StatefulPartitionedCall�
#dense_96656/StatefulPartitionedCallStatefulPartitionedCallxdense_96656_115899084dense_96656_115899086*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96656_layer_call_and_return_conditional_losses_115898986�
1batch_normalization_64528/StatefulPartitionedCallStatefulPartitionedCall,dense_96656/StatefulPartitionedCall:output:0#batch_normalization_64528_115899089#batch_normalization_64528_115899091#batch_normalization_64528_115899093#batch_normalization_64528_115899095*
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
GPU 2J 8� *a
f\RZ
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115898958�
#dense_96657/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_64528/StatefulPartitionedCall:output:0dense_96657_115899098dense_96657_115899100*
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
GPU 2J 8� *S
fNRL
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899012{
IdentityIdentity,dense_96657/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_64528/StatefulPartitionedCall$^dense_96656/StatefulPartitionedCall$^dense_96657/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:����������: : : : : : : : 2f
1batch_normalization_64528/StatefulPartitionedCall1batch_normalization_64528/StatefulPartitionedCall2J
#dense_96656/StatefulPartitionedCall#dense_96656/StatefulPartitionedCall2J
#dense_96657/StatefulPartitionedCall#dense_96657/StatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�
�
/__inference_dense_96656_layer_call_fn_115899350

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
GPU 2J 8� *S
fNRL
J__inference_dense_96656_layer_call_and_return_conditional_losses_115898986o
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
�
�
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899406

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
�	
�
>__inference_shallow_sparse_model_10528_layer_call_fn_115899038
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
GPU 2J 8� *b
f]R[
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899019o
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
_user_specified_name	input_1"�	L
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
>__inference_shallow_sparse_model_10528_layer_call_fn_115899038
>__inference_shallow_sparse_model_10528_layer_call_fn_115899240
>__inference_shallow_sparse_model_10528_layer_call_fn_115899261
>__inference_shallow_sparse_model_10528_layer_call_fn_115899144�
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
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899294
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899341
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899167
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899190�
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
$__inference__wrapped_model_115898887input_1"�
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
@:>	�Z2-shallow_sparse_model_10528/dense_96656/kernel
9:7Z2+shallow_sparse_model_10528/dense_96656/bias
H:FZ2:shallow_sparse_model_10528/batch_normalization_64528/gamma
G:EZ29shallow_sparse_model_10528/batch_normalization_64528/beta
P:NZ (2@shallow_sparse_model_10528/batch_normalization_64528/moving_mean
T:RZ (2Dshallow_sparse_model_10528/batch_normalization_64528/moving_variance
?:=Z2-shallow_sparse_model_10528/dense_96657/kernel
9:72+shallow_sparse_model_10528/dense_96657/bias
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
>__inference_shallow_sparse_model_10528_layer_call_fn_115899038input_1"�
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
>__inference_shallow_sparse_model_10528_layer_call_fn_115899240x"�
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
>__inference_shallow_sparse_model_10528_layer_call_fn_115899261x"�
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
>__inference_shallow_sparse_model_10528_layer_call_fn_115899144input_1"�
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
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899294x"�
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
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899341x"�
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
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899167input_1"�
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
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899190input_1"�
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
/__inference_dense_96656_layer_call_fn_115899350�
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
J__inference_dense_96656_layer_call_and_return_conditional_losses_115899360�
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
=__inference_batch_normalization_64528_layer_call_fn_115899373
=__inference_batch_normalization_64528_layer_call_fn_115899386�
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
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899406
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899440�
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
/__inference_dense_96657_layer_call_fn_115899449�
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
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899460�
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
'__inference_signature_wrapper_115899219input_1"�
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
/__inference_dense_96656_layer_call_fn_115899350inputs"�
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
J__inference_dense_96656_layer_call_and_return_conditional_losses_115899360inputs"�
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
=__inference_batch_normalization_64528_layer_call_fn_115899373inputs"�
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
=__inference_batch_normalization_64528_layer_call_fn_115899386inputs"�
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
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899406inputs"�
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
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899440inputs"�
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
/__inference_dense_96657_layer_call_fn_115899449inputs"�
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
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899460inputs"�
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
E:C	�Z24Adam/shallow_sparse_model_10528/dense_96656/kernel/m
>:<Z22Adam/shallow_sparse_model_10528/dense_96656/bias/m
M:KZ2AAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/m
L:JZ2@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/m
D:BZ24Adam/shallow_sparse_model_10528/dense_96657/kernel/m
>:<22Adam/shallow_sparse_model_10528/dense_96657/bias/m
E:C	�Z24Adam/shallow_sparse_model_10528/dense_96656/kernel/v
>:<Z22Adam/shallow_sparse_model_10528/dense_96656/bias/v
M:KZ2AAdam/shallow_sparse_model_10528/batch_normalization_64528/gamma/v
L:JZ2@Adam/shallow_sparse_model_10528/batch_normalization_64528/beta/v
D:BZ24Adam/shallow_sparse_model_10528/dense_96657/kernel/v
>:<22Adam/shallow_sparse_model_10528/dense_96657/bias/v�
$__inference__wrapped_model_115898887r1�.
'�$
"�
input_1����������
� "3�0
.
output_1"�
output_1����������
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899406b3�0
)�&
 �
inputs���������Z
p 
� "%�"
�
0���������Z
� �
X__inference_batch_normalization_64528_layer_call_and_return_conditional_losses_115899440b3�0
)�&
 �
inputs���������Z
p
� "%�"
�
0���������Z
� �
=__inference_batch_normalization_64528_layer_call_fn_115899373U3�0
)�&
 �
inputs���������Z
p 
� "����������Z�
=__inference_batch_normalization_64528_layer_call_fn_115899386U3�0
)�&
 �
inputs���������Z
p
� "����������Z�
J__inference_dense_96656_layer_call_and_return_conditional_losses_115899360]0�-
&�#
!�
inputs����������
� "%�"
�
0���������Z
� �
/__inference_dense_96656_layer_call_fn_115899350P0�-
&�#
!�
inputs����������
� "����������Z�
J__inference_dense_96657_layer_call_and_return_conditional_losses_115899460\/�,
%�"
 �
inputs���������Z
� "%�"
�
0���������
� �
/__inference_dense_96657_layer_call_fn_115899449O/�,
%�"
 �
inputs���������Z
� "�����������
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899167tA�>
'�$
"�
input_1����������
�

trainingp "%�"
�
0���������
� �
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899190tA�>
'�$
"�
input_1����������
�

trainingp"%�"
�
0���������
� �
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899294n;�8
!�
�
x����������
�

trainingp "%�"
�
0���������
� �
Y__inference_shallow_sparse_model_10528_layer_call_and_return_conditional_losses_115899341n;�8
!�
�
x����������
�

trainingp"%�"
�
0���������
� �
>__inference_shallow_sparse_model_10528_layer_call_fn_115899038gA�>
'�$
"�
input_1����������
�

trainingp "�����������
>__inference_shallow_sparse_model_10528_layer_call_fn_115899144gA�>
'�$
"�
input_1����������
�

trainingp"�����������
>__inference_shallow_sparse_model_10528_layer_call_fn_115899240a;�8
!�
�
x����������
�

trainingp "�����������
>__inference_shallow_sparse_model_10528_layer_call_fn_115899261a;�8
!�
�
x����������
�

trainingp"�����������
'__inference_signature_wrapper_115899219}<�9
� 
2�/
-
input_1"�
input_1����������"3�0
.
output_1"�
output_1���������