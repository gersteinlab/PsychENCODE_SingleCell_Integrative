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
 �"serve*2.10.02v2.10.0-rc3-6-g359c3cdfc5f8��

�
(Adam/dense_model_7064/dense_60074/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*9
shared_name*(Adam/dense_model_7064/dense_60074/bias/v
�
<Adam/dense_model_7064/dense_60074/bias/v/Read/ReadVariableOpReadVariableOp(Adam/dense_model_7064/dense_60074/bias/v*
_output_shapes
:*
dtype0
�
*Adam/dense_model_7064/dense_60074/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*;
shared_name,*Adam/dense_model_7064/dense_60074/kernel/v
�
>Adam/dense_model_7064/dense_60074/kernel/v/Read/ReadVariableOpReadVariableOp*Adam/dense_model_7064/dense_60074/kernel/v*
_output_shapes

:
*
dtype0
�
6Adam/dense_model_7064/batch_normalization_40049/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*G
shared_name86Adam/dense_model_7064/batch_normalization_40049/beta/v
�
JAdam/dense_model_7064/batch_normalization_40049/beta/v/Read/ReadVariableOpReadVariableOp6Adam/dense_model_7064/batch_normalization_40049/beta/v*
_output_shapes
:
*
dtype0
�
7Adam/dense_model_7064/batch_normalization_40049/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*H
shared_name97Adam/dense_model_7064/batch_normalization_40049/gamma/v
�
KAdam/dense_model_7064/batch_normalization_40049/gamma/v/Read/ReadVariableOpReadVariableOp7Adam/dense_model_7064/batch_normalization_40049/gamma/v*
_output_shapes
:
*
dtype0
�
(Adam/dense_model_7064/dense_60073/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*9
shared_name*(Adam/dense_model_7064/dense_60073/bias/v
�
<Adam/dense_model_7064/dense_60073/bias/v/Read/ReadVariableOpReadVariableOp(Adam/dense_model_7064/dense_60073/bias/v*
_output_shapes
:
*
dtype0
�
*Adam/dense_model_7064/dense_60073/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:

*;
shared_name,*Adam/dense_model_7064/dense_60073/kernel/v
�
>Adam/dense_model_7064/dense_60073/kernel/v/Read/ReadVariableOpReadVariableOp*Adam/dense_model_7064/dense_60073/kernel/v*
_output_shapes

:

*
dtype0
�
6Adam/dense_model_7064/batch_normalization_40048/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*G
shared_name86Adam/dense_model_7064/batch_normalization_40048/beta/v
�
JAdam/dense_model_7064/batch_normalization_40048/beta/v/Read/ReadVariableOpReadVariableOp6Adam/dense_model_7064/batch_normalization_40048/beta/v*
_output_shapes
:
*
dtype0
�
7Adam/dense_model_7064/batch_normalization_40048/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*H
shared_name97Adam/dense_model_7064/batch_normalization_40048/gamma/v
�
KAdam/dense_model_7064/batch_normalization_40048/gamma/v/Read/ReadVariableOpReadVariableOp7Adam/dense_model_7064/batch_normalization_40048/gamma/v*
_output_shapes
:
*
dtype0
�
(Adam/dense_model_7064/dense_60072/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*9
shared_name*(Adam/dense_model_7064/dense_60072/bias/v
�
<Adam/dense_model_7064/dense_60072/bias/v/Read/ReadVariableOpReadVariableOp(Adam/dense_model_7064/dense_60072/bias/v*
_output_shapes
:
*
dtype0
�
*Adam/dense_model_7064/dense_60072/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�
*;
shared_name,*Adam/dense_model_7064/dense_60072/kernel/v
�
>Adam/dense_model_7064/dense_60072/kernel/v/Read/ReadVariableOpReadVariableOp*Adam/dense_model_7064/dense_60072/kernel/v*
_output_shapes
:	�
*
dtype0
�
(Adam/dense_model_7064/dense_60074/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*9
shared_name*(Adam/dense_model_7064/dense_60074/bias/m
�
<Adam/dense_model_7064/dense_60074/bias/m/Read/ReadVariableOpReadVariableOp(Adam/dense_model_7064/dense_60074/bias/m*
_output_shapes
:*
dtype0
�
*Adam/dense_model_7064/dense_60074/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*;
shared_name,*Adam/dense_model_7064/dense_60074/kernel/m
�
>Adam/dense_model_7064/dense_60074/kernel/m/Read/ReadVariableOpReadVariableOp*Adam/dense_model_7064/dense_60074/kernel/m*
_output_shapes

:
*
dtype0
�
6Adam/dense_model_7064/batch_normalization_40049/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*G
shared_name86Adam/dense_model_7064/batch_normalization_40049/beta/m
�
JAdam/dense_model_7064/batch_normalization_40049/beta/m/Read/ReadVariableOpReadVariableOp6Adam/dense_model_7064/batch_normalization_40049/beta/m*
_output_shapes
:
*
dtype0
�
7Adam/dense_model_7064/batch_normalization_40049/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*H
shared_name97Adam/dense_model_7064/batch_normalization_40049/gamma/m
�
KAdam/dense_model_7064/batch_normalization_40049/gamma/m/Read/ReadVariableOpReadVariableOp7Adam/dense_model_7064/batch_normalization_40049/gamma/m*
_output_shapes
:
*
dtype0
�
(Adam/dense_model_7064/dense_60073/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*9
shared_name*(Adam/dense_model_7064/dense_60073/bias/m
�
<Adam/dense_model_7064/dense_60073/bias/m/Read/ReadVariableOpReadVariableOp(Adam/dense_model_7064/dense_60073/bias/m*
_output_shapes
:
*
dtype0
�
*Adam/dense_model_7064/dense_60073/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:

*;
shared_name,*Adam/dense_model_7064/dense_60073/kernel/m
�
>Adam/dense_model_7064/dense_60073/kernel/m/Read/ReadVariableOpReadVariableOp*Adam/dense_model_7064/dense_60073/kernel/m*
_output_shapes

:

*
dtype0
�
6Adam/dense_model_7064/batch_normalization_40048/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*G
shared_name86Adam/dense_model_7064/batch_normalization_40048/beta/m
�
JAdam/dense_model_7064/batch_normalization_40048/beta/m/Read/ReadVariableOpReadVariableOp6Adam/dense_model_7064/batch_normalization_40048/beta/m*
_output_shapes
:
*
dtype0
�
7Adam/dense_model_7064/batch_normalization_40048/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*H
shared_name97Adam/dense_model_7064/batch_normalization_40048/gamma/m
�
KAdam/dense_model_7064/batch_normalization_40048/gamma/m/Read/ReadVariableOpReadVariableOp7Adam/dense_model_7064/batch_normalization_40048/gamma/m*
_output_shapes
:
*
dtype0
�
(Adam/dense_model_7064/dense_60072/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*9
shared_name*(Adam/dense_model_7064/dense_60072/bias/m
�
<Adam/dense_model_7064/dense_60072/bias/m/Read/ReadVariableOpReadVariableOp(Adam/dense_model_7064/dense_60072/bias/m*
_output_shapes
:
*
dtype0
�
*Adam/dense_model_7064/dense_60072/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�
*;
shared_name,*Adam/dense_model_7064/dense_60072/kernel/m
�
>Adam/dense_model_7064/dense_60072/kernel/m/Read/ReadVariableOpReadVariableOp*Adam/dense_model_7064/dense_60072/kernel/m*
_output_shapes
:	�
*
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
!dense_model_7064/dense_60074/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!dense_model_7064/dense_60074/bias
�
5dense_model_7064/dense_60074/bias/Read/ReadVariableOpReadVariableOp!dense_model_7064/dense_60074/bias*
_output_shapes
:*
dtype0
�
#dense_model_7064/dense_60074/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:
*4
shared_name%#dense_model_7064/dense_60074/kernel
�
7dense_model_7064/dense_60074/kernel/Read/ReadVariableOpReadVariableOp#dense_model_7064/dense_60074/kernel*
_output_shapes

:
*
dtype0
�
:dense_model_7064/batch_normalization_40049/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*K
shared_name<:dense_model_7064/batch_normalization_40049/moving_variance
�
Ndense_model_7064/batch_normalization_40049/moving_variance/Read/ReadVariableOpReadVariableOp:dense_model_7064/batch_normalization_40049/moving_variance*
_output_shapes
:
*
dtype0
�
6dense_model_7064/batch_normalization_40049/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*G
shared_name86dense_model_7064/batch_normalization_40049/moving_mean
�
Jdense_model_7064/batch_normalization_40049/moving_mean/Read/ReadVariableOpReadVariableOp6dense_model_7064/batch_normalization_40049/moving_mean*
_output_shapes
:
*
dtype0
�
/dense_model_7064/batch_normalization_40049/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*@
shared_name1/dense_model_7064/batch_normalization_40049/beta
�
Cdense_model_7064/batch_normalization_40049/beta/Read/ReadVariableOpReadVariableOp/dense_model_7064/batch_normalization_40049/beta*
_output_shapes
:
*
dtype0
�
0dense_model_7064/batch_normalization_40049/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*A
shared_name20dense_model_7064/batch_normalization_40049/gamma
�
Ddense_model_7064/batch_normalization_40049/gamma/Read/ReadVariableOpReadVariableOp0dense_model_7064/batch_normalization_40049/gamma*
_output_shapes
:
*
dtype0
�
!dense_model_7064/dense_60073/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*2
shared_name#!dense_model_7064/dense_60073/bias
�
5dense_model_7064/dense_60073/bias/Read/ReadVariableOpReadVariableOp!dense_model_7064/dense_60073/bias*
_output_shapes
:
*
dtype0
�
#dense_model_7064/dense_60073/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:

*4
shared_name%#dense_model_7064/dense_60073/kernel
�
7dense_model_7064/dense_60073/kernel/Read/ReadVariableOpReadVariableOp#dense_model_7064/dense_60073/kernel*
_output_shapes

:

*
dtype0
�
:dense_model_7064/batch_normalization_40048/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*K
shared_name<:dense_model_7064/batch_normalization_40048/moving_variance
�
Ndense_model_7064/batch_normalization_40048/moving_variance/Read/ReadVariableOpReadVariableOp:dense_model_7064/batch_normalization_40048/moving_variance*
_output_shapes
:
*
dtype0
�
6dense_model_7064/batch_normalization_40048/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*G
shared_name86dense_model_7064/batch_normalization_40048/moving_mean
�
Jdense_model_7064/batch_normalization_40048/moving_mean/Read/ReadVariableOpReadVariableOp6dense_model_7064/batch_normalization_40048/moving_mean*
_output_shapes
:
*
dtype0
�
/dense_model_7064/batch_normalization_40048/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*@
shared_name1/dense_model_7064/batch_normalization_40048/beta
�
Cdense_model_7064/batch_normalization_40048/beta/Read/ReadVariableOpReadVariableOp/dense_model_7064/batch_normalization_40048/beta*
_output_shapes
:
*
dtype0
�
0dense_model_7064/batch_normalization_40048/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*A
shared_name20dense_model_7064/batch_normalization_40048/gamma
�
Ddense_model_7064/batch_normalization_40048/gamma/Read/ReadVariableOpReadVariableOp0dense_model_7064/batch_normalization_40048/gamma*
_output_shapes
:
*
dtype0
�
!dense_model_7064/dense_60072/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*2
shared_name#!dense_model_7064/dense_60072/bias
�
5dense_model_7064/dense_60072/bias/Read/ReadVariableOpReadVariableOp!dense_model_7064/dense_60072/bias*
_output_shapes
:
*
dtype0
�
#dense_model_7064/dense_60072/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�
*4
shared_name%#dense_model_7064/dense_60072/kernel
�
7dense_model_7064/dense_60072/kernel/Read/ReadVariableOpReadVariableOp#dense_model_7064/dense_60072/kernel*
_output_shapes
:	�
*
dtype0
|
serving_default_input_1Placeholder*(
_output_shapes
:����������*
dtype0*
shape:����������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1#dense_model_7064/dense_60072/kernel!dense_model_7064/dense_60072/bias6dense_model_7064/batch_normalization_40048/moving_mean:dense_model_7064/batch_normalization_40048/moving_variance/dense_model_7064/batch_normalization_40048/beta0dense_model_7064/batch_normalization_40048/gamma#dense_model_7064/dense_60073/kernel!dense_model_7064/dense_60073/bias6dense_model_7064/batch_normalization_40049/moving_mean:dense_model_7064/batch_normalization_40049/moving_variance/dense_model_7064/batch_normalization_40049/beta0dense_model_7064/batch_normalization_40049/gamma#dense_model_7064/dense_60074/kernel!dense_model_7064/dense_60074/bias*
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
GPU 2J 8� */
f*R(
&__inference_signature_wrapper_81179769

NoOpNoOp
�L
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�K
value�KB�K B�K
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
6
"trace_0
#trace_1
$trace_2
%trace_3* 
6
&trace_0
'trace_1
(trace_2
)trace_3* 
* 
�
*	variables
+trainable_variables
,regularization_losses
-	keras_api
.__call__
*/&call_and_return_all_conditional_losses

kernel
bias*
�
0	variables
1trainable_variables
2regularization_losses
3	keras_api
4__call__
*5&call_and_return_all_conditional_losses
6axis
	gamma
beta
moving_mean
moving_variance*
�
7	variables
8trainable_variables
9regularization_losses
:	keras_api
;__call__
*<&call_and_return_all_conditional_losses

kernel
bias*
�
=	variables
>trainable_variables
?regularization_losses
@	keras_api
A__call__
*B&call_and_return_all_conditional_losses
Caxis
	gamma
beta
moving_mean
moving_variance*
�
D	variables
Etrainable_variables
Fregularization_losses
G	keras_api
H__call__
*I&call_and_return_all_conditional_losses

kernel
bias*
�
Jiter

Kbeta_1

Lbeta_2
	Mdecay
Nlearning_ratem�m�m�m�m�m�m�m�m�m�v�v�v�v�v�v�v�v�v�v�*

Oserving_default* 
c]
VARIABLE_VALUE#dense_model_7064/dense_60072/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUE!dense_model_7064/dense_60072/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE0dense_model_7064/batch_normalization_40048/gamma&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
oi
VARIABLE_VALUE/dense_model_7064/batch_normalization_40048/beta&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
vp
VARIABLE_VALUE6dense_model_7064/batch_normalization_40048/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
zt
VARIABLE_VALUE:dense_model_7064/batch_normalization_40048/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
c]
VARIABLE_VALUE#dense_model_7064/dense_60073/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUE!dense_model_7064/dense_60073/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE0dense_model_7064/batch_normalization_40049/gamma&variables/8/.ATTRIBUTES/VARIABLE_VALUE*
oi
VARIABLE_VALUE/dense_model_7064/batch_normalization_40049/beta&variables/9/.ATTRIBUTES/VARIABLE_VALUE*
wq
VARIABLE_VALUE6dense_model_7064/batch_normalization_40049/moving_mean'variables/10/.ATTRIBUTES/VARIABLE_VALUE*
{u
VARIABLE_VALUE:dense_model_7064/batch_normalization_40049/moving_variance'variables/11/.ATTRIBUTES/VARIABLE_VALUE*
d^
VARIABLE_VALUE#dense_model_7064/dense_60074/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE*
b\
VARIABLE_VALUE!dense_model_7064/dense_60074/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE*
 
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
P0
Q1*
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
0
1*

0
1*
* 
�
Rnon_trainable_variables

Slayers
Tmetrics
Ulayer_regularization_losses
Vlayer_metrics
*	variables
+trainable_variables
,regularization_losses
.__call__
*/&call_and_return_all_conditional_losses
&/"call_and_return_conditional_losses*

Wtrace_0* 

Xtrace_0* 
 
0
1
2
3*

0
1*
* 
�
Ynon_trainable_variables

Zlayers
[metrics
\layer_regularization_losses
]layer_metrics
0	variables
1trainable_variables
2regularization_losses
4__call__
*5&call_and_return_all_conditional_losses
&5"call_and_return_conditional_losses*

^trace_0
_trace_1* 

`trace_0
atrace_1* 
* 

0
1*

0
1*
* 
�
bnon_trainable_variables

clayers
dmetrics
elayer_regularization_losses
flayer_metrics
7	variables
8trainable_variables
9regularization_losses
;__call__
*<&call_and_return_all_conditional_losses
&<"call_and_return_conditional_losses*

gtrace_0* 

htrace_0* 
 
0
1
2
3*

0
1*
* 
�
inon_trainable_variables

jlayers
kmetrics
llayer_regularization_losses
mlayer_metrics
=	variables
>trainable_variables
?regularization_losses
A__call__
*B&call_and_return_all_conditional_losses
&B"call_and_return_conditional_losses*

ntrace_0
otrace_1* 

ptrace_0
qtrace_1* 
* 

0
1*

0
1*
* 
�
rnon_trainable_variables

slayers
tmetrics
ulayer_regularization_losses
vlayer_metrics
D	variables
Etrainable_variables
Fregularization_losses
H__call__
*I&call_and_return_all_conditional_losses
&I"call_and_return_conditional_losses*

wtrace_0* 

xtrace_0* 
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
y	variables
z	keras_api
	{total
	|count*
J
}	variables
~	keras_api
	total

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

{0
|1*

y	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

0
�1*

}	variables*
SM
VARIABLE_VALUEtotal4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE*
* 
��
VARIABLE_VALUE*Adam/dense_model_7064/dense_60072/kernel/mBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUE(Adam/dense_model_7064/dense_60072/bias/mBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE7Adam/dense_model_7064/batch_normalization_40048/gamma/mBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE6Adam/dense_model_7064/batch_normalization_40048/beta/mBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE*Adam/dense_model_7064/dense_60073/kernel/mBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUE(Adam/dense_model_7064/dense_60073/bias/mBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE7Adam/dense_model_7064/batch_normalization_40049/gamma/mBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE6Adam/dense_model_7064/batch_normalization_40049/beta/mBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE*Adam/dense_model_7064/dense_60074/kernel/mCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�
VARIABLE_VALUE(Adam/dense_model_7064/dense_60074/bias/mCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE*Adam/dense_model_7064/dense_60072/kernel/vBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUE(Adam/dense_model_7064/dense_60072/bias/vBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE7Adam/dense_model_7064/batch_normalization_40048/gamma/vBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE6Adam/dense_model_7064/batch_normalization_40048/beta/vBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE*Adam/dense_model_7064/dense_60073/kernel/vBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUE(Adam/dense_model_7064/dense_60073/bias/vBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE7Adam/dense_model_7064/batch_normalization_40049/gamma/vBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE6Adam/dense_model_7064/batch_normalization_40049/beta/vBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
��
VARIABLE_VALUE*Adam/dense_model_7064/dense_60074/kernel/vCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�
VARIABLE_VALUE(Adam/dense_model_7064/dense_60074/bias/vCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename7dense_model_7064/dense_60072/kernel/Read/ReadVariableOp5dense_model_7064/dense_60072/bias/Read/ReadVariableOpDdense_model_7064/batch_normalization_40048/gamma/Read/ReadVariableOpCdense_model_7064/batch_normalization_40048/beta/Read/ReadVariableOpJdense_model_7064/batch_normalization_40048/moving_mean/Read/ReadVariableOpNdense_model_7064/batch_normalization_40048/moving_variance/Read/ReadVariableOp7dense_model_7064/dense_60073/kernel/Read/ReadVariableOp5dense_model_7064/dense_60073/bias/Read/ReadVariableOpDdense_model_7064/batch_normalization_40049/gamma/Read/ReadVariableOpCdense_model_7064/batch_normalization_40049/beta/Read/ReadVariableOpJdense_model_7064/batch_normalization_40049/moving_mean/Read/ReadVariableOpNdense_model_7064/batch_normalization_40049/moving_variance/Read/ReadVariableOp7dense_model_7064/dense_60074/kernel/Read/ReadVariableOp5dense_model_7064/dense_60074/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOp>Adam/dense_model_7064/dense_60072/kernel/m/Read/ReadVariableOp<Adam/dense_model_7064/dense_60072/bias/m/Read/ReadVariableOpKAdam/dense_model_7064/batch_normalization_40048/gamma/m/Read/ReadVariableOpJAdam/dense_model_7064/batch_normalization_40048/beta/m/Read/ReadVariableOp>Adam/dense_model_7064/dense_60073/kernel/m/Read/ReadVariableOp<Adam/dense_model_7064/dense_60073/bias/m/Read/ReadVariableOpKAdam/dense_model_7064/batch_normalization_40049/gamma/m/Read/ReadVariableOpJAdam/dense_model_7064/batch_normalization_40049/beta/m/Read/ReadVariableOp>Adam/dense_model_7064/dense_60074/kernel/m/Read/ReadVariableOp<Adam/dense_model_7064/dense_60074/bias/m/Read/ReadVariableOp>Adam/dense_model_7064/dense_60072/kernel/v/Read/ReadVariableOp<Adam/dense_model_7064/dense_60072/bias/v/Read/ReadVariableOpKAdam/dense_model_7064/batch_normalization_40048/gamma/v/Read/ReadVariableOpJAdam/dense_model_7064/batch_normalization_40048/beta/v/Read/ReadVariableOp>Adam/dense_model_7064/dense_60073/kernel/v/Read/ReadVariableOp<Adam/dense_model_7064/dense_60073/bias/v/Read/ReadVariableOpKAdam/dense_model_7064/batch_normalization_40049/gamma/v/Read/ReadVariableOpJAdam/dense_model_7064/batch_normalization_40049/beta/v/Read/ReadVariableOp>Adam/dense_model_7064/dense_60074/kernel/v/Read/ReadVariableOp<Adam/dense_model_7064/dense_60074/bias/v/Read/ReadVariableOpConst*8
Tin1
/2-	*
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
!__inference__traced_save_81180343
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename#dense_model_7064/dense_60072/kernel!dense_model_7064/dense_60072/bias0dense_model_7064/batch_normalization_40048/gamma/dense_model_7064/batch_normalization_40048/beta6dense_model_7064/batch_normalization_40048/moving_mean:dense_model_7064/batch_normalization_40048/moving_variance#dense_model_7064/dense_60073/kernel!dense_model_7064/dense_60073/bias0dense_model_7064/batch_normalization_40049/gamma/dense_model_7064/batch_normalization_40049/beta6dense_model_7064/batch_normalization_40049/moving_mean:dense_model_7064/batch_normalization_40049/moving_variance#dense_model_7064/dense_60074/kernel!dense_model_7064/dense_60074/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_ratetotal_1count_1totalcount*Adam/dense_model_7064/dense_60072/kernel/m(Adam/dense_model_7064/dense_60072/bias/m7Adam/dense_model_7064/batch_normalization_40048/gamma/m6Adam/dense_model_7064/batch_normalization_40048/beta/m*Adam/dense_model_7064/dense_60073/kernel/m(Adam/dense_model_7064/dense_60073/bias/m7Adam/dense_model_7064/batch_normalization_40049/gamma/m6Adam/dense_model_7064/batch_normalization_40049/beta/m*Adam/dense_model_7064/dense_60074/kernel/m(Adam/dense_model_7064/dense_60074/bias/m*Adam/dense_model_7064/dense_60072/kernel/v(Adam/dense_model_7064/dense_60072/bias/v7Adam/dense_model_7064/batch_normalization_40048/gamma/v6Adam/dense_model_7064/batch_normalization_40048/beta/v*Adam/dense_model_7064/dense_60073/kernel/v(Adam/dense_model_7064/dense_60073/bias/v7Adam/dense_model_7064/batch_normalization_40049/gamma/v6Adam/dense_model_7064/batch_normalization_40049/beta/v*Adam/dense_model_7064/dense_60074/kernel/v(Adam/dense_model_7064/dense_60074/bias/v*7
Tin0
.2,*
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
$__inference__traced_restore_81180482��
�$
�
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179371

inputs5
'assignmovingavg_readvariableop_resource:
7
)assignmovingavg_1_readvariableop_resource:
*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:

identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:
�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������
l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:
*
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
:
*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:
x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:
�
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
:
*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:
~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:
�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:
k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�	
�
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179992

inputs1
matmul_readvariableop_resource:	�
-
biasadd_readvariableop_resource:

identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�
*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
w
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
�`
�
!__inference__traced_save_81180343
file_prefixB
>savev2_dense_model_7064_dense_60072_kernel_read_readvariableop@
<savev2_dense_model_7064_dense_60072_bias_read_readvariableopO
Ksavev2_dense_model_7064_batch_normalization_40048_gamma_read_readvariableopN
Jsavev2_dense_model_7064_batch_normalization_40048_beta_read_readvariableopU
Qsavev2_dense_model_7064_batch_normalization_40048_moving_mean_read_readvariableopY
Usavev2_dense_model_7064_batch_normalization_40048_moving_variance_read_readvariableopB
>savev2_dense_model_7064_dense_60073_kernel_read_readvariableop@
<savev2_dense_model_7064_dense_60073_bias_read_readvariableopO
Ksavev2_dense_model_7064_batch_normalization_40049_gamma_read_readvariableopN
Jsavev2_dense_model_7064_batch_normalization_40049_beta_read_readvariableopU
Qsavev2_dense_model_7064_batch_normalization_40049_moving_mean_read_readvariableopY
Usavev2_dense_model_7064_batch_normalization_40049_moving_variance_read_readvariableopB
>savev2_dense_model_7064_dense_60074_kernel_read_readvariableop@
<savev2_dense_model_7064_dense_60074_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableopI
Esavev2_adam_dense_model_7064_dense_60072_kernel_m_read_readvariableopG
Csavev2_adam_dense_model_7064_dense_60072_bias_m_read_readvariableopV
Rsavev2_adam_dense_model_7064_batch_normalization_40048_gamma_m_read_readvariableopU
Qsavev2_adam_dense_model_7064_batch_normalization_40048_beta_m_read_readvariableopI
Esavev2_adam_dense_model_7064_dense_60073_kernel_m_read_readvariableopG
Csavev2_adam_dense_model_7064_dense_60073_bias_m_read_readvariableopV
Rsavev2_adam_dense_model_7064_batch_normalization_40049_gamma_m_read_readvariableopU
Qsavev2_adam_dense_model_7064_batch_normalization_40049_beta_m_read_readvariableopI
Esavev2_adam_dense_model_7064_dense_60074_kernel_m_read_readvariableopG
Csavev2_adam_dense_model_7064_dense_60074_bias_m_read_readvariableopI
Esavev2_adam_dense_model_7064_dense_60072_kernel_v_read_readvariableopG
Csavev2_adam_dense_model_7064_dense_60072_bias_v_read_readvariableopV
Rsavev2_adam_dense_model_7064_batch_normalization_40048_gamma_v_read_readvariableopU
Qsavev2_adam_dense_model_7064_batch_normalization_40048_beta_v_read_readvariableopI
Esavev2_adam_dense_model_7064_dense_60073_kernel_v_read_readvariableopG
Csavev2_adam_dense_model_7064_dense_60073_bias_v_read_readvariableopV
Rsavev2_adam_dense_model_7064_batch_normalization_40049_gamma_v_read_readvariableopU
Qsavev2_adam_dense_model_7064_batch_normalization_40049_beta_v_read_readvariableopI
Esavev2_adam_dense_model_7064_dense_60074_kernel_v_read_readvariableopG
Csavev2_adam_dense_model_7064_dense_60074_bias_v_read_readvariableop
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
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:,*
dtype0*�
value�B�,B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:,*
dtype0*k
valuebB`,B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0>savev2_dense_model_7064_dense_60072_kernel_read_readvariableop<savev2_dense_model_7064_dense_60072_bias_read_readvariableopKsavev2_dense_model_7064_batch_normalization_40048_gamma_read_readvariableopJsavev2_dense_model_7064_batch_normalization_40048_beta_read_readvariableopQsavev2_dense_model_7064_batch_normalization_40048_moving_mean_read_readvariableopUsavev2_dense_model_7064_batch_normalization_40048_moving_variance_read_readvariableop>savev2_dense_model_7064_dense_60073_kernel_read_readvariableop<savev2_dense_model_7064_dense_60073_bias_read_readvariableopKsavev2_dense_model_7064_batch_normalization_40049_gamma_read_readvariableopJsavev2_dense_model_7064_batch_normalization_40049_beta_read_readvariableopQsavev2_dense_model_7064_batch_normalization_40049_moving_mean_read_readvariableopUsavev2_dense_model_7064_batch_normalization_40049_moving_variance_read_readvariableop>savev2_dense_model_7064_dense_60074_kernel_read_readvariableop<savev2_dense_model_7064_dense_60074_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopEsavev2_adam_dense_model_7064_dense_60072_kernel_m_read_readvariableopCsavev2_adam_dense_model_7064_dense_60072_bias_m_read_readvariableopRsavev2_adam_dense_model_7064_batch_normalization_40048_gamma_m_read_readvariableopQsavev2_adam_dense_model_7064_batch_normalization_40048_beta_m_read_readvariableopEsavev2_adam_dense_model_7064_dense_60073_kernel_m_read_readvariableopCsavev2_adam_dense_model_7064_dense_60073_bias_m_read_readvariableopRsavev2_adam_dense_model_7064_batch_normalization_40049_gamma_m_read_readvariableopQsavev2_adam_dense_model_7064_batch_normalization_40049_beta_m_read_readvariableopEsavev2_adam_dense_model_7064_dense_60074_kernel_m_read_readvariableopCsavev2_adam_dense_model_7064_dense_60074_bias_m_read_readvariableopEsavev2_adam_dense_model_7064_dense_60072_kernel_v_read_readvariableopCsavev2_adam_dense_model_7064_dense_60072_bias_v_read_readvariableopRsavev2_adam_dense_model_7064_batch_normalization_40048_gamma_v_read_readvariableopQsavev2_adam_dense_model_7064_batch_normalization_40048_beta_v_read_readvariableopEsavev2_adam_dense_model_7064_dense_60073_kernel_v_read_readvariableopCsavev2_adam_dense_model_7064_dense_60073_bias_v_read_readvariableopRsavev2_adam_dense_model_7064_batch_normalization_40049_gamma_v_read_readvariableopQsavev2_adam_dense_model_7064_batch_normalization_40049_beta_v_read_readvariableopEsavev2_adam_dense_model_7064_dense_60074_kernel_v_read_readvariableopCsavev2_adam_dense_model_7064_dense_60074_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *:
dtypes0
.2,	�
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

identity_1Identity_1:output:0*�
_input_shapes�
�: :	�
:
:
:
:
:
:

:
:
:
:
:
:
:: : : : : : : : : :	�
:
:
:
:

:
:
:
:
::	�
:
:
:
:

:
:
:
:
:: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:%!

_output_shapes
:	�
: 

_output_shapes
:
: 

_output_shapes
:
: 

_output_shapes
:
: 

_output_shapes
:
: 

_output_shapes
:
:$ 

_output_shapes

:

: 

_output_shapes
:
: 	

_output_shapes
:
: 


_output_shapes
:
: 

_output_shapes
:
: 

_output_shapes
:
:$ 

_output_shapes

:
: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :%!

_output_shapes
:	�
: 

_output_shapes
:
: 

_output_shapes
:
: 

_output_shapes
:
:$ 

_output_shapes

:

: 

_output_shapes
:
: 

_output_shapes
:
: 

_output_shapes
:
:$  

_output_shapes

:
: !

_output_shapes
::%"!

_output_shapes
:	�
: #

_output_shapes
:
: $

_output_shapes
:
: %

_output_shapes
:
:$& 

_output_shapes

:

: '

_output_shapes
:
: (

_output_shapes
:
: )

_output_shapes
:
:$* 

_output_shapes

:
: +

_output_shapes
::,

_output_shapes
: 
�
�
<__inference_batch_normalization_40049_layer_call_fn_81180117

inputs
unknown:

	unknown_0:

	unknown_1:

	unknown_2:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179371o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
&__inference_signature_wrapper_81179769
input_1
unknown:	�

	unknown_0:

	unknown_1:

	unknown_2:

	unknown_3:

	unknown_4:

	unknown_5:


	unknown_6:

	unknown_7:

	unknown_8:

	unknown_9:


unknown_10:


unknown_11:


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
GPU 2J 8� *,
f'R%
#__inference__wrapped_model_81179218o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
.__inference_dense_60074_layer_call_fn_81180180

inputs
unknown:

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
I__inference_dense_60074_layer_call_and_return_conditional_losses_81179450o
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
:���������
: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179324

inputs*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:
,
cast_2_readvariableop_resource:
,
cast_3_readvariableop_resource:

identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:
m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�"
�
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179590
x'
dense_60072_81179556:	�
"
dense_60072_81179558:
0
"batch_normalization_40048_81179561:
0
"batch_normalization_40048_81179563:
0
"batch_normalization_40048_81179565:
0
"batch_normalization_40048_81179567:
&
dense_60073_81179570:

"
dense_60073_81179572:
0
"batch_normalization_40049_81179575:
0
"batch_normalization_40049_81179577:
0
"batch_normalization_40049_81179579:
0
"batch_normalization_40049_81179581:
&
dense_60074_81179584:
"
dense_60074_81179586:
identity��1batch_normalization_40048/StatefulPartitionedCall�1batch_normalization_40049/StatefulPartitionedCall�#dense_60072/StatefulPartitionedCall�#dense_60073/StatefulPartitionedCall�#dense_60074/StatefulPartitionedCall�
#dense_60072/StatefulPartitionedCallStatefulPartitionedCallxdense_60072_81179556dense_60072_81179558*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179399�
1batch_normalization_40048/StatefulPartitionedCallStatefulPartitionedCall,dense_60072/StatefulPartitionedCall:output:0"batch_normalization_40048_81179561"batch_normalization_40048_81179563"batch_normalization_40048_81179565"batch_normalization_40048_81179567*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179289�
#dense_60073/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40048/StatefulPartitionedCall:output:0dense_60073_81179570dense_60073_81179572*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60073_layer_call_and_return_conditional_losses_81179424�
1batch_normalization_40049/StatefulPartitionedCallStatefulPartitionedCall,dense_60073/StatefulPartitionedCall:output:0"batch_normalization_40049_81179575"batch_normalization_40049_81179577"batch_normalization_40049_81179579"batch_normalization_40049_81179581*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179371�
#dense_60074/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40049/StatefulPartitionedCall:output:0dense_60074_81179584dense_60074_81179586*
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
I__inference_dense_60074_layer_call_and_return_conditional_losses_81179450{
IdentityIdentity,dense_60074/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_40048/StatefulPartitionedCall2^batch_normalization_40049/StatefulPartitionedCall$^dense_60072/StatefulPartitionedCall$^dense_60073/StatefulPartitionedCall$^dense_60074/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2f
1batch_normalization_40048/StatefulPartitionedCall1batch_normalization_40048/StatefulPartitionedCall2f
1batch_normalization_40049/StatefulPartitionedCall1batch_normalization_40049/StatefulPartitionedCall2J
#dense_60072/StatefulPartitionedCall#dense_60072/StatefulPartitionedCall2J
#dense_60073/StatefulPartitionedCall#dense_60073/StatefulPartitionedCall2J
#dense_60074/StatefulPartitionedCall#dense_60074/StatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�	
�
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179399

inputs1
matmul_readvariableop_resource:	�
-
biasadd_readvariableop_resource:

identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�
*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
w
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
�
3__inference_dense_model_7064_layer_call_fn_81179488
input_1
unknown:	�

	unknown_0:

	unknown_1:

	unknown_2:

	unknown_3:

	unknown_4:

	unknown_5:


	unknown_6:

	unknown_7:

	unknown_8:

	unknown_9:


unknown_10:


unknown_11:


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
GPU 2J 8� *W
fRRP
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179457o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
3__inference_dense_model_7064_layer_call_fn_81179835
x
unknown:	�

	unknown_0:

	unknown_1:

	unknown_2:

	unknown_3:

	unknown_4:

	unknown_5:


	unknown_6:

	unknown_7:

	unknown_8:

	unknown_9:


unknown_10:


unknown_11:


unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
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
GPU 2J 8� *W
fRRP
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179590o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�
�
<__inference_batch_normalization_40048_layer_call_fn_81180018

inputs
unknown:

	unknown_0:

	unknown_1:

	unknown_2:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179289o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180038

inputs*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:
,
cast_2_readvariableop_resource:
,
cast_3_readvariableop_resource:

identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:
m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180137

inputs*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:
,
cast_2_readvariableop_resource:
,
cast_3_readvariableop_resource:

identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:
m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�M
�
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179890
x=
*dense_60072_matmul_readvariableop_resource:	�
9
+dense_60072_biasadd_readvariableop_resource:
D
6batch_normalization_40048_cast_readvariableop_resource:
F
8batch_normalization_40048_cast_1_readvariableop_resource:
F
8batch_normalization_40048_cast_2_readvariableop_resource:
F
8batch_normalization_40048_cast_3_readvariableop_resource:
<
*dense_60073_matmul_readvariableop_resource:

9
+dense_60073_biasadd_readvariableop_resource:
D
6batch_normalization_40049_cast_readvariableop_resource:
F
8batch_normalization_40049_cast_1_readvariableop_resource:
F
8batch_normalization_40049_cast_2_readvariableop_resource:
F
8batch_normalization_40049_cast_3_readvariableop_resource:
<
*dense_60074_matmul_readvariableop_resource:
9
+dense_60074_biasadd_readvariableop_resource:
identity��-batch_normalization_40048/Cast/ReadVariableOp�/batch_normalization_40048/Cast_1/ReadVariableOp�/batch_normalization_40048/Cast_2/ReadVariableOp�/batch_normalization_40048/Cast_3/ReadVariableOp�-batch_normalization_40049/Cast/ReadVariableOp�/batch_normalization_40049/Cast_1/ReadVariableOp�/batch_normalization_40049/Cast_2/ReadVariableOp�/batch_normalization_40049/Cast_3/ReadVariableOp�"dense_60072/BiasAdd/ReadVariableOp�!dense_60072/MatMul/ReadVariableOp�"dense_60073/BiasAdd/ReadVariableOp�!dense_60073/MatMul/ReadVariableOp�"dense_60074/BiasAdd/ReadVariableOp�!dense_60074/MatMul/ReadVariableOp�
!dense_60072/MatMul/ReadVariableOpReadVariableOp*dense_60072_matmul_readvariableop_resource*
_output_shapes
:	�
*
dtype0|
dense_60072/MatMulMatMulx)dense_60072/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
"dense_60072/BiasAdd/ReadVariableOpReadVariableOp+dense_60072_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
dense_60072/BiasAddBiasAdddense_60072/MatMul:product:0*dense_60072/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
-batch_normalization_40048/Cast/ReadVariableOpReadVariableOp6batch_normalization_40048_cast_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40048/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_40048_cast_1_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40048/Cast_2/ReadVariableOpReadVariableOp8batch_normalization_40048_cast_2_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40048/Cast_3/ReadVariableOpReadVariableOp8batch_normalization_40048_cast_3_readvariableop_resource*
_output_shapes
:
*
dtype0n
)batch_normalization_40048/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_40048/batchnorm/addAddV27batch_normalization_40048/Cast_1/ReadVariableOp:value:02batch_normalization_40048/batchnorm/add/y:output:0*
T0*
_output_shapes
:
�
)batch_normalization_40048/batchnorm/RsqrtRsqrt+batch_normalization_40048/batchnorm/add:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40048/batchnorm/mulMul-batch_normalization_40048/batchnorm/Rsqrt:y:07batch_normalization_40048/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
�
)batch_normalization_40048/batchnorm/mul_1Muldense_60072/BiasAdd:output:0+batch_normalization_40048/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������
�
)batch_normalization_40048/batchnorm/mul_2Mul5batch_normalization_40048/Cast/ReadVariableOp:value:0+batch_normalization_40048/batchnorm/mul:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40048/batchnorm/subSub7batch_normalization_40048/Cast_2/ReadVariableOp:value:0-batch_normalization_40048/batchnorm/mul_2:z:0*
T0*
_output_shapes
:
�
)batch_normalization_40048/batchnorm/add_1AddV2-batch_normalization_40048/batchnorm/mul_1:z:0+batch_normalization_40048/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
�
!dense_60073/MatMul/ReadVariableOpReadVariableOp*dense_60073_matmul_readvariableop_resource*
_output_shapes

:

*
dtype0�
dense_60073/MatMulMatMul-batch_normalization_40048/batchnorm/add_1:z:0)dense_60073/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
"dense_60073/BiasAdd/ReadVariableOpReadVariableOp+dense_60073_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
dense_60073/BiasAddBiasAdddense_60073/MatMul:product:0*dense_60073/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
-batch_normalization_40049/Cast/ReadVariableOpReadVariableOp6batch_normalization_40049_cast_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40049/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_40049_cast_1_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40049/Cast_2/ReadVariableOpReadVariableOp8batch_normalization_40049_cast_2_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40049/Cast_3/ReadVariableOpReadVariableOp8batch_normalization_40049_cast_3_readvariableop_resource*
_output_shapes
:
*
dtype0n
)batch_normalization_40049/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_40049/batchnorm/addAddV27batch_normalization_40049/Cast_1/ReadVariableOp:value:02batch_normalization_40049/batchnorm/add/y:output:0*
T0*
_output_shapes
:
�
)batch_normalization_40049/batchnorm/RsqrtRsqrt+batch_normalization_40049/batchnorm/add:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40049/batchnorm/mulMul-batch_normalization_40049/batchnorm/Rsqrt:y:07batch_normalization_40049/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
�
)batch_normalization_40049/batchnorm/mul_1Muldense_60073/BiasAdd:output:0+batch_normalization_40049/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������
�
)batch_normalization_40049/batchnorm/mul_2Mul5batch_normalization_40049/Cast/ReadVariableOp:value:0+batch_normalization_40049/batchnorm/mul:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40049/batchnorm/subSub7batch_normalization_40049/Cast_2/ReadVariableOp:value:0-batch_normalization_40049/batchnorm/mul_2:z:0*
T0*
_output_shapes
:
�
)batch_normalization_40049/batchnorm/add_1AddV2-batch_normalization_40049/batchnorm/mul_1:z:0+batch_normalization_40049/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
�
!dense_60074/MatMul/ReadVariableOpReadVariableOp*dense_60074_matmul_readvariableop_resource*
_output_shapes

:
*
dtype0�
dense_60074/MatMulMatMul-batch_normalization_40049/batchnorm/add_1:z:0)dense_60074/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
"dense_60074/BiasAdd/ReadVariableOpReadVariableOp+dense_60074_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_60074/BiasAddBiasAdddense_60074/MatMul:product:0*dense_60074/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������n
dense_60074/SigmoidSigmoiddense_60074/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
IdentityIdentitydense_60074/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp.^batch_normalization_40048/Cast/ReadVariableOp0^batch_normalization_40048/Cast_1/ReadVariableOp0^batch_normalization_40048/Cast_2/ReadVariableOp0^batch_normalization_40048/Cast_3/ReadVariableOp.^batch_normalization_40049/Cast/ReadVariableOp0^batch_normalization_40049/Cast_1/ReadVariableOp0^batch_normalization_40049/Cast_2/ReadVariableOp0^batch_normalization_40049/Cast_3/ReadVariableOp#^dense_60072/BiasAdd/ReadVariableOp"^dense_60072/MatMul/ReadVariableOp#^dense_60073/BiasAdd/ReadVariableOp"^dense_60073/MatMul/ReadVariableOp#^dense_60074/BiasAdd/ReadVariableOp"^dense_60074/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2^
-batch_normalization_40048/Cast/ReadVariableOp-batch_normalization_40048/Cast/ReadVariableOp2b
/batch_normalization_40048/Cast_1/ReadVariableOp/batch_normalization_40048/Cast_1/ReadVariableOp2b
/batch_normalization_40048/Cast_2/ReadVariableOp/batch_normalization_40048/Cast_2/ReadVariableOp2b
/batch_normalization_40048/Cast_3/ReadVariableOp/batch_normalization_40048/Cast_3/ReadVariableOp2^
-batch_normalization_40049/Cast/ReadVariableOp-batch_normalization_40049/Cast/ReadVariableOp2b
/batch_normalization_40049/Cast_1/ReadVariableOp/batch_normalization_40049/Cast_1/ReadVariableOp2b
/batch_normalization_40049/Cast_2/ReadVariableOp/batch_normalization_40049/Cast_2/ReadVariableOp2b
/batch_normalization_40049/Cast_3/ReadVariableOp/batch_normalization_40049/Cast_3/ReadVariableOp2H
"dense_60072/BiasAdd/ReadVariableOp"dense_60072/BiasAdd/ReadVariableOp2F
!dense_60072/MatMul/ReadVariableOp!dense_60072/MatMul/ReadVariableOp2H
"dense_60073/BiasAdd/ReadVariableOp"dense_60073/BiasAdd/ReadVariableOp2F
!dense_60073/MatMul/ReadVariableOp!dense_60073/MatMul/ReadVariableOp2H
"dense_60074/BiasAdd/ReadVariableOp"dense_60074/BiasAdd/ReadVariableOp2F
!dense_60074/MatMul/ReadVariableOp!dense_60074/MatMul/ReadVariableOp:K G
(
_output_shapes
:����������

_user_specified_namex
�$
�
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180171

inputs5
'assignmovingavg_readvariableop_resource:
7
)assignmovingavg_1_readvariableop_resource:
*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:

identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:
�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������
l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:
*
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
:
*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:
x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:
�
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
:
*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:
~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:
�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:
k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�$
�
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179289

inputs5
'assignmovingavg_readvariableop_resource:
7
)assignmovingavg_1_readvariableop_resource:
*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:

identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:
�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������
l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:
*
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
:
*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:
x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:
�
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
:
*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:
~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:
�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:
k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
3__inference_dense_model_7064_layer_call_fn_81179654
input_1
unknown:	�

	unknown_0:

	unknown_1:

	unknown_2:

	unknown_3:

	unknown_4:

	unknown_5:


	unknown_6:

	unknown_7:

	unknown_8:

	unknown_9:


unknown_10:


unknown_11:


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
GPU 2J 8� *W
fRRP
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179590o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�"
�
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179728
input_1'
dense_60072_81179694:	�
"
dense_60072_81179696:
0
"batch_normalization_40048_81179699:
0
"batch_normalization_40048_81179701:
0
"batch_normalization_40048_81179703:
0
"batch_normalization_40048_81179705:
&
dense_60073_81179708:

"
dense_60073_81179710:
0
"batch_normalization_40049_81179713:
0
"batch_normalization_40049_81179715:
0
"batch_normalization_40049_81179717:
0
"batch_normalization_40049_81179719:
&
dense_60074_81179722:
"
dense_60074_81179724:
identity��1batch_normalization_40048/StatefulPartitionedCall�1batch_normalization_40049/StatefulPartitionedCall�#dense_60072/StatefulPartitionedCall�#dense_60073/StatefulPartitionedCall�#dense_60074/StatefulPartitionedCall�
#dense_60072/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_60072_81179694dense_60072_81179696*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179399�
1batch_normalization_40048/StatefulPartitionedCallStatefulPartitionedCall,dense_60072/StatefulPartitionedCall:output:0"batch_normalization_40048_81179699"batch_normalization_40048_81179701"batch_normalization_40048_81179703"batch_normalization_40048_81179705*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179289�
#dense_60073/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40048/StatefulPartitionedCall:output:0dense_60073_81179708dense_60073_81179710*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60073_layer_call_and_return_conditional_losses_81179424�
1batch_normalization_40049/StatefulPartitionedCallStatefulPartitionedCall,dense_60073/StatefulPartitionedCall:output:0"batch_normalization_40049_81179713"batch_normalization_40049_81179715"batch_normalization_40049_81179717"batch_normalization_40049_81179719*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179371�
#dense_60074/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40049/StatefulPartitionedCall:output:0dense_60074_81179722dense_60074_81179724*
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
I__inference_dense_60074_layer_call_and_return_conditional_losses_81179450{
IdentityIdentity,dense_60074/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_40048/StatefulPartitionedCall2^batch_normalization_40049/StatefulPartitionedCall$^dense_60072/StatefulPartitionedCall$^dense_60073/StatefulPartitionedCall$^dense_60074/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2f
1batch_normalization_40048/StatefulPartitionedCall1batch_normalization_40048/StatefulPartitionedCall2f
1batch_normalization_40049/StatefulPartitionedCall1batch_normalization_40049/StatefulPartitionedCall2J
#dense_60072/StatefulPartitionedCall#dense_60072/StatefulPartitionedCall2J
#dense_60073/StatefulPartitionedCall#dense_60073/StatefulPartitionedCall2J
#dense_60074/StatefulPartitionedCall#dense_60074/StatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�"
�
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179457
x'
dense_60072_81179400:	�
"
dense_60072_81179402:
0
"batch_normalization_40048_81179405:
0
"batch_normalization_40048_81179407:
0
"batch_normalization_40048_81179409:
0
"batch_normalization_40048_81179411:
&
dense_60073_81179425:

"
dense_60073_81179427:
0
"batch_normalization_40049_81179430:
0
"batch_normalization_40049_81179432:
0
"batch_normalization_40049_81179434:
0
"batch_normalization_40049_81179436:
&
dense_60074_81179451:
"
dense_60074_81179453:
identity��1batch_normalization_40048/StatefulPartitionedCall�1batch_normalization_40049/StatefulPartitionedCall�#dense_60072/StatefulPartitionedCall�#dense_60073/StatefulPartitionedCall�#dense_60074/StatefulPartitionedCall�
#dense_60072/StatefulPartitionedCallStatefulPartitionedCallxdense_60072_81179400dense_60072_81179402*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179399�
1batch_normalization_40048/StatefulPartitionedCallStatefulPartitionedCall,dense_60072/StatefulPartitionedCall:output:0"batch_normalization_40048_81179405"batch_normalization_40048_81179407"batch_normalization_40048_81179409"batch_normalization_40048_81179411*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179242�
#dense_60073/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40048/StatefulPartitionedCall:output:0dense_60073_81179425dense_60073_81179427*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60073_layer_call_and_return_conditional_losses_81179424�
1batch_normalization_40049/StatefulPartitionedCallStatefulPartitionedCall,dense_60073/StatefulPartitionedCall:output:0"batch_normalization_40049_81179430"batch_normalization_40049_81179432"batch_normalization_40049_81179434"batch_normalization_40049_81179436*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179324�
#dense_60074/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40049/StatefulPartitionedCall:output:0dense_60074_81179451dense_60074_81179453*
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
I__inference_dense_60074_layer_call_and_return_conditional_losses_81179450{
IdentityIdentity,dense_60074/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_40048/StatefulPartitionedCall2^batch_normalization_40049/StatefulPartitionedCall$^dense_60072/StatefulPartitionedCall$^dense_60073/StatefulPartitionedCall$^dense_60074/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2f
1batch_normalization_40048/StatefulPartitionedCall1batch_normalization_40048/StatefulPartitionedCall2f
1batch_normalization_40049/StatefulPartitionedCall1batch_normalization_40049/StatefulPartitionedCall2J
#dense_60072/StatefulPartitionedCall#dense_60072/StatefulPartitionedCall2J
#dense_60073/StatefulPartitionedCall#dense_60073/StatefulPartitionedCall2J
#dense_60074/StatefulPartitionedCall#dense_60074/StatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�
�
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179242

inputs*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:
,
cast_2_readvariableop_resource:
,
cast_3_readvariableop_resource:

identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:
m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp2.
Cast_2/ReadVariableOpCast_2/ReadVariableOp2.
Cast_3/ReadVariableOpCast_3/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
3__inference_dense_model_7064_layer_call_fn_81179802
x
unknown:	�

	unknown_0:

	unknown_1:

	unknown_2:

	unknown_3:

	unknown_4:

	unknown_5:


	unknown_6:

	unknown_7:

	unknown_8:

	unknown_9:


unknown_10:


unknown_11:


unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
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
GPU 2J 8� *W
fRRP
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179457o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:K G
(
_output_shapes
:����������

_user_specified_namex
�b
�
#__inference__wrapped_model_81179218
input_1N
;dense_model_7064_dense_60072_matmul_readvariableop_resource:	�
J
<dense_model_7064_dense_60072_biasadd_readvariableop_resource:
U
Gdense_model_7064_batch_normalization_40048_cast_readvariableop_resource:
W
Idense_model_7064_batch_normalization_40048_cast_1_readvariableop_resource:
W
Idense_model_7064_batch_normalization_40048_cast_2_readvariableop_resource:
W
Idense_model_7064_batch_normalization_40048_cast_3_readvariableop_resource:
M
;dense_model_7064_dense_60073_matmul_readvariableop_resource:

J
<dense_model_7064_dense_60073_biasadd_readvariableop_resource:
U
Gdense_model_7064_batch_normalization_40049_cast_readvariableop_resource:
W
Idense_model_7064_batch_normalization_40049_cast_1_readvariableop_resource:
W
Idense_model_7064_batch_normalization_40049_cast_2_readvariableop_resource:
W
Idense_model_7064_batch_normalization_40049_cast_3_readvariableop_resource:
M
;dense_model_7064_dense_60074_matmul_readvariableop_resource:
J
<dense_model_7064_dense_60074_biasadd_readvariableop_resource:
identity��>dense_model_7064/batch_normalization_40048/Cast/ReadVariableOp�@dense_model_7064/batch_normalization_40048/Cast_1/ReadVariableOp�@dense_model_7064/batch_normalization_40048/Cast_2/ReadVariableOp�@dense_model_7064/batch_normalization_40048/Cast_3/ReadVariableOp�>dense_model_7064/batch_normalization_40049/Cast/ReadVariableOp�@dense_model_7064/batch_normalization_40049/Cast_1/ReadVariableOp�@dense_model_7064/batch_normalization_40049/Cast_2/ReadVariableOp�@dense_model_7064/batch_normalization_40049/Cast_3/ReadVariableOp�3dense_model_7064/dense_60072/BiasAdd/ReadVariableOp�2dense_model_7064/dense_60072/MatMul/ReadVariableOp�3dense_model_7064/dense_60073/BiasAdd/ReadVariableOp�2dense_model_7064/dense_60073/MatMul/ReadVariableOp�3dense_model_7064/dense_60074/BiasAdd/ReadVariableOp�2dense_model_7064/dense_60074/MatMul/ReadVariableOp�
2dense_model_7064/dense_60072/MatMul/ReadVariableOpReadVariableOp;dense_model_7064_dense_60072_matmul_readvariableop_resource*
_output_shapes
:	�
*
dtype0�
#dense_model_7064/dense_60072/MatMulMatMulinput_1:dense_model_7064/dense_60072/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
3dense_model_7064/dense_60072/BiasAdd/ReadVariableOpReadVariableOp<dense_model_7064_dense_60072_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
$dense_model_7064/dense_60072/BiasAddBiasAdd-dense_model_7064/dense_60072/MatMul:product:0;dense_model_7064/dense_60072/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
>dense_model_7064/batch_normalization_40048/Cast/ReadVariableOpReadVariableOpGdense_model_7064_batch_normalization_40048_cast_readvariableop_resource*
_output_shapes
:
*
dtype0�
@dense_model_7064/batch_normalization_40048/Cast_1/ReadVariableOpReadVariableOpIdense_model_7064_batch_normalization_40048_cast_1_readvariableop_resource*
_output_shapes
:
*
dtype0�
@dense_model_7064/batch_normalization_40048/Cast_2/ReadVariableOpReadVariableOpIdense_model_7064_batch_normalization_40048_cast_2_readvariableop_resource*
_output_shapes
:
*
dtype0�
@dense_model_7064/batch_normalization_40048/Cast_3/ReadVariableOpReadVariableOpIdense_model_7064_batch_normalization_40048_cast_3_readvariableop_resource*
_output_shapes
:
*
dtype0
:dense_model_7064/batch_normalization_40048/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
8dense_model_7064/batch_normalization_40048/batchnorm/addAddV2Hdense_model_7064/batch_normalization_40048/Cast_1/ReadVariableOp:value:0Cdense_model_7064/batch_normalization_40048/batchnorm/add/y:output:0*
T0*
_output_shapes
:
�
:dense_model_7064/batch_normalization_40048/batchnorm/RsqrtRsqrt<dense_model_7064/batch_normalization_40048/batchnorm/add:z:0*
T0*
_output_shapes
:
�
8dense_model_7064/batch_normalization_40048/batchnorm/mulMul>dense_model_7064/batch_normalization_40048/batchnorm/Rsqrt:y:0Hdense_model_7064/batch_normalization_40048/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
�
:dense_model_7064/batch_normalization_40048/batchnorm/mul_1Mul-dense_model_7064/dense_60072/BiasAdd:output:0<dense_model_7064/batch_normalization_40048/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������
�
:dense_model_7064/batch_normalization_40048/batchnorm/mul_2MulFdense_model_7064/batch_normalization_40048/Cast/ReadVariableOp:value:0<dense_model_7064/batch_normalization_40048/batchnorm/mul:z:0*
T0*
_output_shapes
:
�
8dense_model_7064/batch_normalization_40048/batchnorm/subSubHdense_model_7064/batch_normalization_40048/Cast_2/ReadVariableOp:value:0>dense_model_7064/batch_normalization_40048/batchnorm/mul_2:z:0*
T0*
_output_shapes
:
�
:dense_model_7064/batch_normalization_40048/batchnorm/add_1AddV2>dense_model_7064/batch_normalization_40048/batchnorm/mul_1:z:0<dense_model_7064/batch_normalization_40048/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
�
2dense_model_7064/dense_60073/MatMul/ReadVariableOpReadVariableOp;dense_model_7064_dense_60073_matmul_readvariableop_resource*
_output_shapes

:

*
dtype0�
#dense_model_7064/dense_60073/MatMulMatMul>dense_model_7064/batch_normalization_40048/batchnorm/add_1:z:0:dense_model_7064/dense_60073/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
3dense_model_7064/dense_60073/BiasAdd/ReadVariableOpReadVariableOp<dense_model_7064_dense_60073_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
$dense_model_7064/dense_60073/BiasAddBiasAdd-dense_model_7064/dense_60073/MatMul:product:0;dense_model_7064/dense_60073/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
>dense_model_7064/batch_normalization_40049/Cast/ReadVariableOpReadVariableOpGdense_model_7064_batch_normalization_40049_cast_readvariableop_resource*
_output_shapes
:
*
dtype0�
@dense_model_7064/batch_normalization_40049/Cast_1/ReadVariableOpReadVariableOpIdense_model_7064_batch_normalization_40049_cast_1_readvariableop_resource*
_output_shapes
:
*
dtype0�
@dense_model_7064/batch_normalization_40049/Cast_2/ReadVariableOpReadVariableOpIdense_model_7064_batch_normalization_40049_cast_2_readvariableop_resource*
_output_shapes
:
*
dtype0�
@dense_model_7064/batch_normalization_40049/Cast_3/ReadVariableOpReadVariableOpIdense_model_7064_batch_normalization_40049_cast_3_readvariableop_resource*
_output_shapes
:
*
dtype0
:dense_model_7064/batch_normalization_40049/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
8dense_model_7064/batch_normalization_40049/batchnorm/addAddV2Hdense_model_7064/batch_normalization_40049/Cast_1/ReadVariableOp:value:0Cdense_model_7064/batch_normalization_40049/batchnorm/add/y:output:0*
T0*
_output_shapes
:
�
:dense_model_7064/batch_normalization_40049/batchnorm/RsqrtRsqrt<dense_model_7064/batch_normalization_40049/batchnorm/add:z:0*
T0*
_output_shapes
:
�
8dense_model_7064/batch_normalization_40049/batchnorm/mulMul>dense_model_7064/batch_normalization_40049/batchnorm/Rsqrt:y:0Hdense_model_7064/batch_normalization_40049/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:
�
:dense_model_7064/batch_normalization_40049/batchnorm/mul_1Mul-dense_model_7064/dense_60073/BiasAdd:output:0<dense_model_7064/batch_normalization_40049/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������
�
:dense_model_7064/batch_normalization_40049/batchnorm/mul_2MulFdense_model_7064/batch_normalization_40049/Cast/ReadVariableOp:value:0<dense_model_7064/batch_normalization_40049/batchnorm/mul:z:0*
T0*
_output_shapes
:
�
8dense_model_7064/batch_normalization_40049/batchnorm/subSubHdense_model_7064/batch_normalization_40049/Cast_2/ReadVariableOp:value:0>dense_model_7064/batch_normalization_40049/batchnorm/mul_2:z:0*
T0*
_output_shapes
:
�
:dense_model_7064/batch_normalization_40049/batchnorm/add_1AddV2>dense_model_7064/batch_normalization_40049/batchnorm/mul_1:z:0<dense_model_7064/batch_normalization_40049/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
�
2dense_model_7064/dense_60074/MatMul/ReadVariableOpReadVariableOp;dense_model_7064_dense_60074_matmul_readvariableop_resource*
_output_shapes

:
*
dtype0�
#dense_model_7064/dense_60074/MatMulMatMul>dense_model_7064/batch_normalization_40049/batchnorm/add_1:z:0:dense_model_7064/dense_60074/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
3dense_model_7064/dense_60074/BiasAdd/ReadVariableOpReadVariableOp<dense_model_7064_dense_60074_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
$dense_model_7064/dense_60074/BiasAddBiasAdd-dense_model_7064/dense_60074/MatMul:product:0;dense_model_7064/dense_60074/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
$dense_model_7064/dense_60074/SigmoidSigmoid-dense_model_7064/dense_60074/BiasAdd:output:0*
T0*'
_output_shapes
:���������w
IdentityIdentity(dense_model_7064/dense_60074/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp?^dense_model_7064/batch_normalization_40048/Cast/ReadVariableOpA^dense_model_7064/batch_normalization_40048/Cast_1/ReadVariableOpA^dense_model_7064/batch_normalization_40048/Cast_2/ReadVariableOpA^dense_model_7064/batch_normalization_40048/Cast_3/ReadVariableOp?^dense_model_7064/batch_normalization_40049/Cast/ReadVariableOpA^dense_model_7064/batch_normalization_40049/Cast_1/ReadVariableOpA^dense_model_7064/batch_normalization_40049/Cast_2/ReadVariableOpA^dense_model_7064/batch_normalization_40049/Cast_3/ReadVariableOp4^dense_model_7064/dense_60072/BiasAdd/ReadVariableOp3^dense_model_7064/dense_60072/MatMul/ReadVariableOp4^dense_model_7064/dense_60073/BiasAdd/ReadVariableOp3^dense_model_7064/dense_60073/MatMul/ReadVariableOp4^dense_model_7064/dense_60074/BiasAdd/ReadVariableOp3^dense_model_7064/dense_60074/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2�
>dense_model_7064/batch_normalization_40048/Cast/ReadVariableOp>dense_model_7064/batch_normalization_40048/Cast/ReadVariableOp2�
@dense_model_7064/batch_normalization_40048/Cast_1/ReadVariableOp@dense_model_7064/batch_normalization_40048/Cast_1/ReadVariableOp2�
@dense_model_7064/batch_normalization_40048/Cast_2/ReadVariableOp@dense_model_7064/batch_normalization_40048/Cast_2/ReadVariableOp2�
@dense_model_7064/batch_normalization_40048/Cast_3/ReadVariableOp@dense_model_7064/batch_normalization_40048/Cast_3/ReadVariableOp2�
>dense_model_7064/batch_normalization_40049/Cast/ReadVariableOp>dense_model_7064/batch_normalization_40049/Cast/ReadVariableOp2�
@dense_model_7064/batch_normalization_40049/Cast_1/ReadVariableOp@dense_model_7064/batch_normalization_40049/Cast_1/ReadVariableOp2�
@dense_model_7064/batch_normalization_40049/Cast_2/ReadVariableOp@dense_model_7064/batch_normalization_40049/Cast_2/ReadVariableOp2�
@dense_model_7064/batch_normalization_40049/Cast_3/ReadVariableOp@dense_model_7064/batch_normalization_40049/Cast_3/ReadVariableOp2j
3dense_model_7064/dense_60072/BiasAdd/ReadVariableOp3dense_model_7064/dense_60072/BiasAdd/ReadVariableOp2h
2dense_model_7064/dense_60072/MatMul/ReadVariableOp2dense_model_7064/dense_60072/MatMul/ReadVariableOp2j
3dense_model_7064/dense_60073/BiasAdd/ReadVariableOp3dense_model_7064/dense_60073/BiasAdd/ReadVariableOp2h
2dense_model_7064/dense_60073/MatMul/ReadVariableOp2dense_model_7064/dense_60073/MatMul/ReadVariableOp2j
3dense_model_7064/dense_60074/BiasAdd/ReadVariableOp3dense_model_7064/dense_60074/BiasAdd/ReadVariableOp2h
2dense_model_7064/dense_60074/MatMul/ReadVariableOp2dense_model_7064/dense_60074/MatMul/ReadVariableOp:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
޴
� 
$__inference__traced_restore_81180482
file_prefixG
4assignvariableop_dense_model_7064_dense_60072_kernel:	�
B
4assignvariableop_1_dense_model_7064_dense_60072_bias:
Q
Cassignvariableop_2_dense_model_7064_batch_normalization_40048_gamma:
P
Bassignvariableop_3_dense_model_7064_batch_normalization_40048_beta:
W
Iassignvariableop_4_dense_model_7064_batch_normalization_40048_moving_mean:
[
Massignvariableop_5_dense_model_7064_batch_normalization_40048_moving_variance:
H
6assignvariableop_6_dense_model_7064_dense_60073_kernel:

B
4assignvariableop_7_dense_model_7064_dense_60073_bias:
Q
Cassignvariableop_8_dense_model_7064_batch_normalization_40049_gamma:
P
Bassignvariableop_9_dense_model_7064_batch_normalization_40049_beta:
X
Jassignvariableop_10_dense_model_7064_batch_normalization_40049_moving_mean:
\
Nassignvariableop_11_dense_model_7064_batch_normalization_40049_moving_variance:
I
7assignvariableop_12_dense_model_7064_dense_60074_kernel:
C
5assignvariableop_13_dense_model_7064_dense_60074_bias:'
assignvariableop_14_adam_iter:	 )
assignvariableop_15_adam_beta_1: )
assignvariableop_16_adam_beta_2: (
assignvariableop_17_adam_decay: 0
&assignvariableop_18_adam_learning_rate: %
assignvariableop_19_total_1: %
assignvariableop_20_count_1: #
assignvariableop_21_total: #
assignvariableop_22_count: Q
>assignvariableop_23_adam_dense_model_7064_dense_60072_kernel_m:	�
J
<assignvariableop_24_adam_dense_model_7064_dense_60072_bias_m:
Y
Kassignvariableop_25_adam_dense_model_7064_batch_normalization_40048_gamma_m:
X
Jassignvariableop_26_adam_dense_model_7064_batch_normalization_40048_beta_m:
P
>assignvariableop_27_adam_dense_model_7064_dense_60073_kernel_m:

J
<assignvariableop_28_adam_dense_model_7064_dense_60073_bias_m:
Y
Kassignvariableop_29_adam_dense_model_7064_batch_normalization_40049_gamma_m:
X
Jassignvariableop_30_adam_dense_model_7064_batch_normalization_40049_beta_m:
P
>assignvariableop_31_adam_dense_model_7064_dense_60074_kernel_m:
J
<assignvariableop_32_adam_dense_model_7064_dense_60074_bias_m:Q
>assignvariableop_33_adam_dense_model_7064_dense_60072_kernel_v:	�
J
<assignvariableop_34_adam_dense_model_7064_dense_60072_bias_v:
Y
Kassignvariableop_35_adam_dense_model_7064_batch_normalization_40048_gamma_v:
X
Jassignvariableop_36_adam_dense_model_7064_batch_normalization_40048_beta_v:
P
>assignvariableop_37_adam_dense_model_7064_dense_60073_kernel_v:

J
<assignvariableop_38_adam_dense_model_7064_dense_60073_bias_v:
Y
Kassignvariableop_39_adam_dense_model_7064_batch_normalization_40049_gamma_v:
X
Jassignvariableop_40_adam_dense_model_7064_batch_normalization_40049_beta_v:
P
>assignvariableop_41_adam_dense_model_7064_dense_60074_kernel_v:
J
<assignvariableop_42_adam_dense_model_7064_dense_60074_bias_v:
identity_44��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:,*
dtype0*�
value�B�,B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:,*
dtype0*k
valuebB`,B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�::::::::::::::::::::::::::::::::::::::::::::*:
dtypes0
.2,	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOp4assignvariableop_dense_model_7064_dense_60072_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp4assignvariableop_1_dense_model_7064_dense_60072_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOpCassignvariableop_2_dense_model_7064_batch_normalization_40048_gammaIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpBassignvariableop_3_dense_model_7064_batch_normalization_40048_betaIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpIassignvariableop_4_dense_model_7064_batch_normalization_40048_moving_meanIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpMassignvariableop_5_dense_model_7064_batch_normalization_40048_moving_varianceIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp6assignvariableop_6_dense_model_7064_dense_60073_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp4assignvariableop_7_dense_model_7064_dense_60073_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOpCassignvariableop_8_dense_model_7064_batch_normalization_40049_gammaIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOpBassignvariableop_9_dense_model_7064_batch_normalization_40049_betaIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOpJassignvariableop_10_dense_model_7064_batch_normalization_40049_moving_meanIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpNassignvariableop_11_dense_model_7064_batch_normalization_40049_moving_varianceIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp7assignvariableop_12_dense_model_7064_dense_60074_kernelIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp5assignvariableop_13_dense_model_7064_dense_60074_biasIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_14AssignVariableOpassignvariableop_14_adam_iterIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOpassignvariableop_15_adam_beta_1Identity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOpassignvariableop_16_adam_beta_2Identity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOpassignvariableop_17_adam_decayIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp&assignvariableop_18_adam_learning_rateIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOpassignvariableop_19_total_1Identity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOpassignvariableop_20_count_1Identity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOpassignvariableop_21_totalIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOpassignvariableop_22_countIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOp>assignvariableop_23_adam_dense_model_7064_dense_60072_kernel_mIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOp<assignvariableop_24_adam_dense_model_7064_dense_60072_bias_mIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOpKassignvariableop_25_adam_dense_model_7064_batch_normalization_40048_gamma_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpJassignvariableop_26_adam_dense_model_7064_batch_normalization_40048_beta_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOp>assignvariableop_27_adam_dense_model_7064_dense_60073_kernel_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOp<assignvariableop_28_adam_dense_model_7064_dense_60073_bias_mIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpKassignvariableop_29_adam_dense_model_7064_batch_normalization_40049_gamma_mIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOpJassignvariableop_30_adam_dense_model_7064_batch_normalization_40049_beta_mIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOp>assignvariableop_31_adam_dense_model_7064_dense_60074_kernel_mIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOp<assignvariableop_32_adam_dense_model_7064_dense_60074_bias_mIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOp>assignvariableop_33_adam_dense_model_7064_dense_60072_kernel_vIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOp<assignvariableop_34_adam_dense_model_7064_dense_60072_bias_vIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOpKassignvariableop_35_adam_dense_model_7064_batch_normalization_40048_gamma_vIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_36AssignVariableOpJassignvariableop_36_adam_dense_model_7064_batch_normalization_40048_beta_vIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_37AssignVariableOp>assignvariableop_37_adam_dense_model_7064_dense_60073_kernel_vIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_38AssignVariableOp<assignvariableop_38_adam_dense_model_7064_dense_60073_bias_vIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_39AssignVariableOpKassignvariableop_39_adam_dense_model_7064_batch_normalization_40049_gamma_vIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_40AssignVariableOpJassignvariableop_40_adam_dense_model_7064_batch_normalization_40049_beta_vIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_41AssignVariableOp>assignvariableop_41_adam_dense_model_7064_dense_60074_kernel_vIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_42AssignVariableOp<assignvariableop_42_adam_dense_model_7064_dense_60074_bias_vIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype01
NoOpNoOp"/device:CPU:0*
_output_shapes
 �
Identity_43Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_44IdentityIdentity_43:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_44Identity_44:output:0*k
_input_shapesZ
X: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
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
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
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
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422(
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
�$
�
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180072

inputs5
'assignmovingavg_readvariableop_resource:
7
)assignmovingavg_1_readvariableop_resource:
*
cast_readvariableop_resource:
,
cast_1_readvariableop_resource:

identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:
�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������
l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:
*
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
:
*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:
x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:
�
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
:
*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:
~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:
�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:
*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:
*
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
:
P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:
m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:
c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������
h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:
k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:
r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2*
Cast/ReadVariableOpCast/ReadVariableOp2.
Cast_1/ReadVariableOpCast_1/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
<__inference_batch_normalization_40048_layer_call_fn_81180005

inputs
unknown:

	unknown_0:

	unknown_1:

	unknown_2:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179242o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
ل
�
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179973
x=
*dense_60072_matmul_readvariableop_resource:	�
9
+dense_60072_biasadd_readvariableop_resource:
O
Abatch_normalization_40048_assignmovingavg_readvariableop_resource:
Q
Cbatch_normalization_40048_assignmovingavg_1_readvariableop_resource:
D
6batch_normalization_40048_cast_readvariableop_resource:
F
8batch_normalization_40048_cast_1_readvariableop_resource:
<
*dense_60073_matmul_readvariableop_resource:

9
+dense_60073_biasadd_readvariableop_resource:
O
Abatch_normalization_40049_assignmovingavg_readvariableop_resource:
Q
Cbatch_normalization_40049_assignmovingavg_1_readvariableop_resource:
D
6batch_normalization_40049_cast_readvariableop_resource:
F
8batch_normalization_40049_cast_1_readvariableop_resource:
<
*dense_60074_matmul_readvariableop_resource:
9
+dense_60074_biasadd_readvariableop_resource:
identity��)batch_normalization_40048/AssignMovingAvg�8batch_normalization_40048/AssignMovingAvg/ReadVariableOp�+batch_normalization_40048/AssignMovingAvg_1�:batch_normalization_40048/AssignMovingAvg_1/ReadVariableOp�-batch_normalization_40048/Cast/ReadVariableOp�/batch_normalization_40048/Cast_1/ReadVariableOp�)batch_normalization_40049/AssignMovingAvg�8batch_normalization_40049/AssignMovingAvg/ReadVariableOp�+batch_normalization_40049/AssignMovingAvg_1�:batch_normalization_40049/AssignMovingAvg_1/ReadVariableOp�-batch_normalization_40049/Cast/ReadVariableOp�/batch_normalization_40049/Cast_1/ReadVariableOp�"dense_60072/BiasAdd/ReadVariableOp�!dense_60072/MatMul/ReadVariableOp�"dense_60073/BiasAdd/ReadVariableOp�!dense_60073/MatMul/ReadVariableOp�"dense_60074/BiasAdd/ReadVariableOp�!dense_60074/MatMul/ReadVariableOp�
!dense_60072/MatMul/ReadVariableOpReadVariableOp*dense_60072_matmul_readvariableop_resource*
_output_shapes
:	�
*
dtype0|
dense_60072/MatMulMatMulx)dense_60072/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
"dense_60072/BiasAdd/ReadVariableOpReadVariableOp+dense_60072_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
dense_60072/BiasAddBiasAdddense_60072/MatMul:product:0*dense_60072/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
8batch_normalization_40048/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
&batch_normalization_40048/moments/meanMeandense_60072/BiasAdd:output:0Abatch_normalization_40048/moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(�
.batch_normalization_40048/moments/StopGradientStopGradient/batch_normalization_40048/moments/mean:output:0*
T0*
_output_shapes

:
�
3batch_normalization_40048/moments/SquaredDifferenceSquaredDifferencedense_60072/BiasAdd:output:07batch_normalization_40048/moments/StopGradient:output:0*
T0*'
_output_shapes
:���������
�
<batch_normalization_40048/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
*batch_normalization_40048/moments/varianceMean7batch_normalization_40048/moments/SquaredDifference:z:0Ebatch_normalization_40048/moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(�
)batch_normalization_40048/moments/SqueezeSqueeze/batch_normalization_40048/moments/mean:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 �
+batch_normalization_40048/moments/Squeeze_1Squeeze3batch_normalization_40048/moments/variance:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 t
/batch_normalization_40048/AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
8batch_normalization_40048/AssignMovingAvg/ReadVariableOpReadVariableOpAbatch_normalization_40048_assignmovingavg_readvariableop_resource*
_output_shapes
:
*
dtype0�
-batch_normalization_40048/AssignMovingAvg/subSub@batch_normalization_40048/AssignMovingAvg/ReadVariableOp:value:02batch_normalization_40048/moments/Squeeze:output:0*
T0*
_output_shapes
:
�
-batch_normalization_40048/AssignMovingAvg/mulMul1batch_normalization_40048/AssignMovingAvg/sub:z:08batch_normalization_40048/AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:
�
)batch_normalization_40048/AssignMovingAvgAssignSubVariableOpAbatch_normalization_40048_assignmovingavg_readvariableop_resource1batch_normalization_40048/AssignMovingAvg/mul:z:09^batch_normalization_40048/AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0v
1batch_normalization_40048/AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
:batch_normalization_40048/AssignMovingAvg_1/ReadVariableOpReadVariableOpCbatch_normalization_40048_assignmovingavg_1_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40048/AssignMovingAvg_1/subSubBbatch_normalization_40048/AssignMovingAvg_1/ReadVariableOp:value:04batch_normalization_40048/moments/Squeeze_1:output:0*
T0*
_output_shapes
:
�
/batch_normalization_40048/AssignMovingAvg_1/mulMul3batch_normalization_40048/AssignMovingAvg_1/sub:z:0:batch_normalization_40048/AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:
�
+batch_normalization_40048/AssignMovingAvg_1AssignSubVariableOpCbatch_normalization_40048_assignmovingavg_1_readvariableop_resource3batch_normalization_40048/AssignMovingAvg_1/mul:z:0;^batch_normalization_40048/AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0�
-batch_normalization_40048/Cast/ReadVariableOpReadVariableOp6batch_normalization_40048_cast_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40048/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_40048_cast_1_readvariableop_resource*
_output_shapes
:
*
dtype0n
)batch_normalization_40048/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_40048/batchnorm/addAddV24batch_normalization_40048/moments/Squeeze_1:output:02batch_normalization_40048/batchnorm/add/y:output:0*
T0*
_output_shapes
:
�
)batch_normalization_40048/batchnorm/RsqrtRsqrt+batch_normalization_40048/batchnorm/add:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40048/batchnorm/mulMul-batch_normalization_40048/batchnorm/Rsqrt:y:07batch_normalization_40048/Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:
�
)batch_normalization_40048/batchnorm/mul_1Muldense_60072/BiasAdd:output:0+batch_normalization_40048/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������
�
)batch_normalization_40048/batchnorm/mul_2Mul2batch_normalization_40048/moments/Squeeze:output:0+batch_normalization_40048/batchnorm/mul:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40048/batchnorm/subSub5batch_normalization_40048/Cast/ReadVariableOp:value:0-batch_normalization_40048/batchnorm/mul_2:z:0*
T0*
_output_shapes
:
�
)batch_normalization_40048/batchnorm/add_1AddV2-batch_normalization_40048/batchnorm/mul_1:z:0+batch_normalization_40048/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
�
!dense_60073/MatMul/ReadVariableOpReadVariableOp*dense_60073_matmul_readvariableop_resource*
_output_shapes

:

*
dtype0�
dense_60073/MatMulMatMul-batch_normalization_40048/batchnorm/add_1:z:0)dense_60073/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
"dense_60073/BiasAdd/ReadVariableOpReadVariableOp+dense_60073_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
dense_60073/BiasAddBiasAdddense_60073/MatMul:product:0*dense_60073/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
8batch_normalization_40049/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
&batch_normalization_40049/moments/meanMeandense_60073/BiasAdd:output:0Abatch_normalization_40049/moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(�
.batch_normalization_40049/moments/StopGradientStopGradient/batch_normalization_40049/moments/mean:output:0*
T0*
_output_shapes

:
�
3batch_normalization_40049/moments/SquaredDifferenceSquaredDifferencedense_60073/BiasAdd:output:07batch_normalization_40049/moments/StopGradient:output:0*
T0*'
_output_shapes
:���������
�
<batch_normalization_40049/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
*batch_normalization_40049/moments/varianceMean7batch_normalization_40049/moments/SquaredDifference:z:0Ebatch_normalization_40049/moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:
*
	keep_dims(�
)batch_normalization_40049/moments/SqueezeSqueeze/batch_normalization_40049/moments/mean:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 �
+batch_normalization_40049/moments/Squeeze_1Squeeze3batch_normalization_40049/moments/variance:output:0*
T0*
_output_shapes
:
*
squeeze_dims
 t
/batch_normalization_40049/AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
8batch_normalization_40049/AssignMovingAvg/ReadVariableOpReadVariableOpAbatch_normalization_40049_assignmovingavg_readvariableop_resource*
_output_shapes
:
*
dtype0�
-batch_normalization_40049/AssignMovingAvg/subSub@batch_normalization_40049/AssignMovingAvg/ReadVariableOp:value:02batch_normalization_40049/moments/Squeeze:output:0*
T0*
_output_shapes
:
�
-batch_normalization_40049/AssignMovingAvg/mulMul1batch_normalization_40049/AssignMovingAvg/sub:z:08batch_normalization_40049/AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:
�
)batch_normalization_40049/AssignMovingAvgAssignSubVariableOpAbatch_normalization_40049_assignmovingavg_readvariableop_resource1batch_normalization_40049/AssignMovingAvg/mul:z:09^batch_normalization_40049/AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype0v
1batch_normalization_40049/AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
:batch_normalization_40049/AssignMovingAvg_1/ReadVariableOpReadVariableOpCbatch_normalization_40049_assignmovingavg_1_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40049/AssignMovingAvg_1/subSubBbatch_normalization_40049/AssignMovingAvg_1/ReadVariableOp:value:04batch_normalization_40049/moments/Squeeze_1:output:0*
T0*
_output_shapes
:
�
/batch_normalization_40049/AssignMovingAvg_1/mulMul3batch_normalization_40049/AssignMovingAvg_1/sub:z:0:batch_normalization_40049/AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:
�
+batch_normalization_40049/AssignMovingAvg_1AssignSubVariableOpCbatch_normalization_40049_assignmovingavg_1_readvariableop_resource3batch_normalization_40049/AssignMovingAvg_1/mul:z:0;^batch_normalization_40049/AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0�
-batch_normalization_40049/Cast/ReadVariableOpReadVariableOp6batch_normalization_40049_cast_readvariableop_resource*
_output_shapes
:
*
dtype0�
/batch_normalization_40049/Cast_1/ReadVariableOpReadVariableOp8batch_normalization_40049_cast_1_readvariableop_resource*
_output_shapes
:
*
dtype0n
)batch_normalization_40049/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
'batch_normalization_40049/batchnorm/addAddV24batch_normalization_40049/moments/Squeeze_1:output:02batch_normalization_40049/batchnorm/add/y:output:0*
T0*
_output_shapes
:
�
)batch_normalization_40049/batchnorm/RsqrtRsqrt+batch_normalization_40049/batchnorm/add:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40049/batchnorm/mulMul-batch_normalization_40049/batchnorm/Rsqrt:y:07batch_normalization_40049/Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:
�
)batch_normalization_40049/batchnorm/mul_1Muldense_60073/BiasAdd:output:0+batch_normalization_40049/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������
�
)batch_normalization_40049/batchnorm/mul_2Mul2batch_normalization_40049/moments/Squeeze:output:0+batch_normalization_40049/batchnorm/mul:z:0*
T0*
_output_shapes
:
�
'batch_normalization_40049/batchnorm/subSub5batch_normalization_40049/Cast/ReadVariableOp:value:0-batch_normalization_40049/batchnorm/mul_2:z:0*
T0*
_output_shapes
:
�
)batch_normalization_40049/batchnorm/add_1AddV2-batch_normalization_40049/batchnorm/mul_1:z:0+batch_normalization_40049/batchnorm/sub:z:0*
T0*'
_output_shapes
:���������
�
!dense_60074/MatMul/ReadVariableOpReadVariableOp*dense_60074_matmul_readvariableop_resource*
_output_shapes

:
*
dtype0�
dense_60074/MatMulMatMul-batch_normalization_40049/batchnorm/add_1:z:0)dense_60074/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
"dense_60074/BiasAdd/ReadVariableOpReadVariableOp+dense_60074_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_60074/BiasAddBiasAdddense_60074/MatMul:product:0*dense_60074/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������n
dense_60074/SigmoidSigmoiddense_60074/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
IdentityIdentitydense_60074/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp*^batch_normalization_40048/AssignMovingAvg9^batch_normalization_40048/AssignMovingAvg/ReadVariableOp,^batch_normalization_40048/AssignMovingAvg_1;^batch_normalization_40048/AssignMovingAvg_1/ReadVariableOp.^batch_normalization_40048/Cast/ReadVariableOp0^batch_normalization_40048/Cast_1/ReadVariableOp*^batch_normalization_40049/AssignMovingAvg9^batch_normalization_40049/AssignMovingAvg/ReadVariableOp,^batch_normalization_40049/AssignMovingAvg_1;^batch_normalization_40049/AssignMovingAvg_1/ReadVariableOp.^batch_normalization_40049/Cast/ReadVariableOp0^batch_normalization_40049/Cast_1/ReadVariableOp#^dense_60072/BiasAdd/ReadVariableOp"^dense_60072/MatMul/ReadVariableOp#^dense_60073/BiasAdd/ReadVariableOp"^dense_60073/MatMul/ReadVariableOp#^dense_60074/BiasAdd/ReadVariableOp"^dense_60074/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2V
)batch_normalization_40048/AssignMovingAvg)batch_normalization_40048/AssignMovingAvg2t
8batch_normalization_40048/AssignMovingAvg/ReadVariableOp8batch_normalization_40048/AssignMovingAvg/ReadVariableOp2Z
+batch_normalization_40048/AssignMovingAvg_1+batch_normalization_40048/AssignMovingAvg_12x
:batch_normalization_40048/AssignMovingAvg_1/ReadVariableOp:batch_normalization_40048/AssignMovingAvg_1/ReadVariableOp2^
-batch_normalization_40048/Cast/ReadVariableOp-batch_normalization_40048/Cast/ReadVariableOp2b
/batch_normalization_40048/Cast_1/ReadVariableOp/batch_normalization_40048/Cast_1/ReadVariableOp2V
)batch_normalization_40049/AssignMovingAvg)batch_normalization_40049/AssignMovingAvg2t
8batch_normalization_40049/AssignMovingAvg/ReadVariableOp8batch_normalization_40049/AssignMovingAvg/ReadVariableOp2Z
+batch_normalization_40049/AssignMovingAvg_1+batch_normalization_40049/AssignMovingAvg_12x
:batch_normalization_40049/AssignMovingAvg_1/ReadVariableOp:batch_normalization_40049/AssignMovingAvg_1/ReadVariableOp2^
-batch_normalization_40049/Cast/ReadVariableOp-batch_normalization_40049/Cast/ReadVariableOp2b
/batch_normalization_40049/Cast_1/ReadVariableOp/batch_normalization_40049/Cast_1/ReadVariableOp2H
"dense_60072/BiasAdd/ReadVariableOp"dense_60072/BiasAdd/ReadVariableOp2F
!dense_60072/MatMul/ReadVariableOp!dense_60072/MatMul/ReadVariableOp2H
"dense_60073/BiasAdd/ReadVariableOp"dense_60073/BiasAdd/ReadVariableOp2F
!dense_60073/MatMul/ReadVariableOp!dense_60073/MatMul/ReadVariableOp2H
"dense_60074/BiasAdd/ReadVariableOp"dense_60074/BiasAdd/ReadVariableOp2F
!dense_60074/MatMul/ReadVariableOp!dense_60074/MatMul/ReadVariableOp:K G
(
_output_shapes
:����������

_user_specified_namex
�	
�
I__inference_dense_60073_layer_call_and_return_conditional_losses_81179424

inputs0
matmul_readvariableop_resource:

-
biasadd_readvariableop_resource:

identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:

*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������
: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�

�
I__inference_dense_60074_layer_call_and_return_conditional_losses_81180191

inputs0
matmul_readvariableop_resource:
-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
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
:���������
: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�

�
I__inference_dense_60074_layer_call_and_return_conditional_losses_81179450

inputs0
matmul_readvariableop_resource:
-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:
*
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
:���������
: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�	
�
I__inference_dense_60073_layer_call_and_return_conditional_losses_81180091

inputs0
matmul_readvariableop_resource:

-
biasadd_readvariableop_resource:

identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:

*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������
: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
<__inference_batch_normalization_40049_layer_call_fn_81180104

inputs
unknown:

	unknown_0:

	unknown_1:

	unknown_2:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179324o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�
�
.__inference_dense_60073_layer_call_fn_81180081

inputs
unknown:


	unknown_0:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60073_layer_call_and_return_conditional_losses_81179424o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������
: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������

 
_user_specified_nameinputs
�"
�
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179691
input_1'
dense_60072_81179657:	�
"
dense_60072_81179659:
0
"batch_normalization_40048_81179662:
0
"batch_normalization_40048_81179664:
0
"batch_normalization_40048_81179666:
0
"batch_normalization_40048_81179668:
&
dense_60073_81179671:

"
dense_60073_81179673:
0
"batch_normalization_40049_81179676:
0
"batch_normalization_40049_81179678:
0
"batch_normalization_40049_81179680:
0
"batch_normalization_40049_81179682:
&
dense_60074_81179685:
"
dense_60074_81179687:
identity��1batch_normalization_40048/StatefulPartitionedCall�1batch_normalization_40049/StatefulPartitionedCall�#dense_60072/StatefulPartitionedCall�#dense_60073/StatefulPartitionedCall�#dense_60074/StatefulPartitionedCall�
#dense_60072/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_60072_81179657dense_60072_81179659*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179399�
1batch_normalization_40048/StatefulPartitionedCallStatefulPartitionedCall,dense_60072/StatefulPartitionedCall:output:0"batch_normalization_40048_81179662"batch_normalization_40048_81179664"batch_normalization_40048_81179666"batch_normalization_40048_81179668*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81179242�
#dense_60073/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40048/StatefulPartitionedCall:output:0dense_60073_81179671dense_60073_81179673*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60073_layer_call_and_return_conditional_losses_81179424�
1batch_normalization_40049/StatefulPartitionedCallStatefulPartitionedCall,dense_60073/StatefulPartitionedCall:output:0"batch_normalization_40049_81179676"batch_normalization_40049_81179678"batch_normalization_40049_81179680"batch_normalization_40049_81179682*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *`
f[RY
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81179324�
#dense_60074/StatefulPartitionedCallStatefulPartitionedCall:batch_normalization_40049/StatefulPartitionedCall:output:0dense_60074_81179685dense_60074_81179687*
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
I__inference_dense_60074_layer_call_and_return_conditional_losses_81179450{
IdentityIdentity,dense_60074/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp2^batch_normalization_40048/StatefulPartitionedCall2^batch_normalization_40049/StatefulPartitionedCall$^dense_60072/StatefulPartitionedCall$^dense_60073/StatefulPartitionedCall$^dense_60074/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*C
_input_shapes2
0:����������: : : : : : : : : : : : : : 2f
1batch_normalization_40048/StatefulPartitionedCall1batch_normalization_40048/StatefulPartitionedCall2f
1batch_normalization_40049/StatefulPartitionedCall1batch_normalization_40049/StatefulPartitionedCall2J
#dense_60072/StatefulPartitionedCall#dense_60072/StatefulPartitionedCall2J
#dense_60073/StatefulPartitionedCall#dense_60073/StatefulPartitionedCall2J
#dense_60074/StatefulPartitionedCall#dense_60074/StatefulPartitionedCall:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
.__inference_dense_60072_layer_call_fn_81179982

inputs
unknown:	�

	unknown_0:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179399o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
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
�
"trace_0
#trace_1
$trace_2
%trace_32�
3__inference_dense_model_7064_layer_call_fn_81179488
3__inference_dense_model_7064_layer_call_fn_81179802
3__inference_dense_model_7064_layer_call_fn_81179835
3__inference_dense_model_7064_layer_call_fn_81179654�
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
 z"trace_0z#trace_1z$trace_2z%trace_3
�
&trace_0
'trace_1
(trace_2
)trace_32�
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179890
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179973
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179691
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179728�
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
 z&trace_0z'trace_1z(trace_2z)trace_3
�B�
#__inference__wrapped_model_81179218input_1"�
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
*	variables
+trainable_variables
,regularization_losses
-	keras_api
.__call__
*/&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
0	variables
1trainable_variables
2regularization_losses
3	keras_api
4__call__
*5&call_and_return_all_conditional_losses
6axis
	gamma
beta
moving_mean
moving_variance"
_tf_keras_layer
�
7	variables
8trainable_variables
9regularization_losses
:	keras_api
;__call__
*<&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
=	variables
>trainable_variables
?regularization_losses
@	keras_api
A__call__
*B&call_and_return_all_conditional_losses
Caxis
	gamma
beta
moving_mean
moving_variance"
_tf_keras_layer
�
D	variables
Etrainable_variables
Fregularization_losses
G	keras_api
H__call__
*I&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
Jiter

Kbeta_1

Lbeta_2
	Mdecay
Nlearning_ratem�m�m�m�m�m�m�m�m�m�v�v�v�v�v�v�v�v�v�v�"
	optimizer
,
Oserving_default"
signature_map
6:4	�
2#dense_model_7064/dense_60072/kernel
/:-
2!dense_model_7064/dense_60072/bias
>:<
20dense_model_7064/batch_normalization_40048/gamma
=:;
2/dense_model_7064/batch_normalization_40048/beta
F:D
 (26dense_model_7064/batch_normalization_40048/moving_mean
J:H
 (2:dense_model_7064/batch_normalization_40048/moving_variance
5:3

2#dense_model_7064/dense_60073/kernel
/:-
2!dense_model_7064/dense_60073/bias
>:<
20dense_model_7064/batch_normalization_40049/gamma
=:;
2/dense_model_7064/batch_normalization_40049/beta
F:D
 (26dense_model_7064/batch_normalization_40049/moving_mean
J:H
 (2:dense_model_7064/batch_normalization_40049/moving_variance
5:3
2#dense_model_7064/dense_60074/kernel
/:-2!dense_model_7064/dense_60074/bias
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
P0
Q1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
3__inference_dense_model_7064_layer_call_fn_81179488input_1"�
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
�B�
3__inference_dense_model_7064_layer_call_fn_81179802x"�
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
�B�
3__inference_dense_model_7064_layer_call_fn_81179835x"�
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
�B�
3__inference_dense_model_7064_layer_call_fn_81179654input_1"�
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
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179890x"�
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
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179973x"�
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
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179691input_1"�
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
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179728input_1"�
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
Rnon_trainable_variables

Slayers
Tmetrics
Ulayer_regularization_losses
Vlayer_metrics
*	variables
+trainable_variables
,regularization_losses
.__call__
*/&call_and_return_all_conditional_losses
&/"call_and_return_conditional_losses"
_generic_user_object
�
Wtrace_02�
.__inference_dense_60072_layer_call_fn_81179982�
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
 zWtrace_0
�
Xtrace_02�
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179992�
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
 zXtrace_0
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
Ynon_trainable_variables

Zlayers
[metrics
\layer_regularization_losses
]layer_metrics
0	variables
1trainable_variables
2regularization_losses
4__call__
*5&call_and_return_all_conditional_losses
&5"call_and_return_conditional_losses"
_generic_user_object
�
^trace_0
_trace_12�
<__inference_batch_normalization_40048_layer_call_fn_81180005
<__inference_batch_normalization_40048_layer_call_fn_81180018�
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
 z^trace_0z_trace_1
�
`trace_0
atrace_12�
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180038
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180072�
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
 z`trace_0zatrace_1
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
bnon_trainable_variables

clayers
dmetrics
elayer_regularization_losses
flayer_metrics
7	variables
8trainable_variables
9regularization_losses
;__call__
*<&call_and_return_all_conditional_losses
&<"call_and_return_conditional_losses"
_generic_user_object
�
gtrace_02�
.__inference_dense_60073_layer_call_fn_81180081�
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
 zgtrace_0
�
htrace_02�
I__inference_dense_60073_layer_call_and_return_conditional_losses_81180091�
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
 zhtrace_0
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
inon_trainable_variables

jlayers
kmetrics
llayer_regularization_losses
mlayer_metrics
=	variables
>trainable_variables
?regularization_losses
A__call__
*B&call_and_return_all_conditional_losses
&B"call_and_return_conditional_losses"
_generic_user_object
�
ntrace_0
otrace_12�
<__inference_batch_normalization_40049_layer_call_fn_81180104
<__inference_batch_normalization_40049_layer_call_fn_81180117�
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
 zntrace_0zotrace_1
�
ptrace_0
qtrace_12�
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180137
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180171�
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
 zptrace_0zqtrace_1
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
rnon_trainable_variables

slayers
tmetrics
ulayer_regularization_losses
vlayer_metrics
D	variables
Etrainable_variables
Fregularization_losses
H__call__
*I&call_and_return_all_conditional_losses
&I"call_and_return_conditional_losses"
_generic_user_object
�
wtrace_02�
.__inference_dense_60074_layer_call_fn_81180180�
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
 zwtrace_0
�
xtrace_02�
I__inference_dense_60074_layer_call_and_return_conditional_losses_81180191�
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
 zxtrace_0
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
�B�
&__inference_signature_wrapper_81179769input_1"�
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
y	variables
z	keras_api
	{total
	|count"
_tf_keras_metric
`
}	variables
~	keras_api
	total

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
.__inference_dense_60072_layer_call_fn_81179982inputs"�
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
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179992inputs"�
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
�B�
<__inference_batch_normalization_40048_layer_call_fn_81180005inputs"�
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
<__inference_batch_normalization_40048_layer_call_fn_81180018inputs"�
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
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180038inputs"�
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
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180072inputs"�
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
.__inference_dense_60073_layer_call_fn_81180081inputs"�
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
I__inference_dense_60073_layer_call_and_return_conditional_losses_81180091inputs"�
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
�B�
<__inference_batch_normalization_40049_layer_call_fn_81180104inputs"�
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
<__inference_batch_normalization_40049_layer_call_fn_81180117inputs"�
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
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180137inputs"�
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
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180171inputs"�
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
.__inference_dense_60074_layer_call_fn_81180180inputs"�
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
I__inference_dense_60074_layer_call_and_return_conditional_losses_81180191inputs"�
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
{0
|1"
trackable_list_wrapper
-
y	variables"
_generic_user_object
:  (2total
:  (2count
/
0
�1"
trackable_list_wrapper
-
}	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
;:9	�
2*Adam/dense_model_7064/dense_60072/kernel/m
4:2
2(Adam/dense_model_7064/dense_60072/bias/m
C:A
27Adam/dense_model_7064/batch_normalization_40048/gamma/m
B:@
26Adam/dense_model_7064/batch_normalization_40048/beta/m
::8

2*Adam/dense_model_7064/dense_60073/kernel/m
4:2
2(Adam/dense_model_7064/dense_60073/bias/m
C:A
27Adam/dense_model_7064/batch_normalization_40049/gamma/m
B:@
26Adam/dense_model_7064/batch_normalization_40049/beta/m
::8
2*Adam/dense_model_7064/dense_60074/kernel/m
4:22(Adam/dense_model_7064/dense_60074/bias/m
;:9	�
2*Adam/dense_model_7064/dense_60072/kernel/v
4:2
2(Adam/dense_model_7064/dense_60072/bias/v
C:A
27Adam/dense_model_7064/batch_normalization_40048/gamma/v
B:@
26Adam/dense_model_7064/batch_normalization_40048/beta/v
::8

2*Adam/dense_model_7064/dense_60073/kernel/v
4:2
2(Adam/dense_model_7064/dense_60073/bias/v
C:A
27Adam/dense_model_7064/batch_normalization_40049/gamma/v
B:@
26Adam/dense_model_7064/batch_normalization_40049/beta/v
::8
2*Adam/dense_model_7064/dense_60074/kernel/v
4:22(Adam/dense_model_7064/dense_60074/bias/v�
#__inference__wrapped_model_81179218x1�.
'�$
"�
input_1����������
� "3�0
.
output_1"�
output_1����������
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180038b3�0
)�&
 �
inputs���������

p 
� "%�"
�
0���������

� �
W__inference_batch_normalization_40048_layer_call_and_return_conditional_losses_81180072b3�0
)�&
 �
inputs���������

p
� "%�"
�
0���������

� �
<__inference_batch_normalization_40048_layer_call_fn_81180005U3�0
)�&
 �
inputs���������

p 
� "����������
�
<__inference_batch_normalization_40048_layer_call_fn_81180018U3�0
)�&
 �
inputs���������

p
� "����������
�
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180137b3�0
)�&
 �
inputs���������

p 
� "%�"
�
0���������

� �
W__inference_batch_normalization_40049_layer_call_and_return_conditional_losses_81180171b3�0
)�&
 �
inputs���������

p
� "%�"
�
0���������

� �
<__inference_batch_normalization_40049_layer_call_fn_81180104U3�0
)�&
 �
inputs���������

p 
� "����������
�
<__inference_batch_normalization_40049_layer_call_fn_81180117U3�0
)�&
 �
inputs���������

p
� "����������
�
I__inference_dense_60072_layer_call_and_return_conditional_losses_81179992]0�-
&�#
!�
inputs����������
� "%�"
�
0���������

� �
.__inference_dense_60072_layer_call_fn_81179982P0�-
&�#
!�
inputs����������
� "����������
�
I__inference_dense_60073_layer_call_and_return_conditional_losses_81180091\/�,
%�"
 �
inputs���������

� "%�"
�
0���������

� �
.__inference_dense_60073_layer_call_fn_81180081O/�,
%�"
 �
inputs���������

� "����������
�
I__inference_dense_60074_layer_call_and_return_conditional_losses_81180191\/�,
%�"
 �
inputs���������

� "%�"
�
0���������
� �
.__inference_dense_60074_layer_call_fn_81180180O/�,
%�"
 �
inputs���������

� "�����������
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179691zA�>
'�$
"�
input_1����������
�

trainingp "%�"
�
0���������
� �
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179728zA�>
'�$
"�
input_1����������
�

trainingp"%�"
�
0���������
� �
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179890t;�8
!�
�
x����������
�

trainingp "%�"
�
0���������
� �
N__inference_dense_model_7064_layer_call_and_return_conditional_losses_81179973t;�8
!�
�
x����������
�

trainingp"%�"
�
0���������
� �
3__inference_dense_model_7064_layer_call_fn_81179488mA�>
'�$
"�
input_1����������
�

trainingp "�����������
3__inference_dense_model_7064_layer_call_fn_81179654mA�>
'�$
"�
input_1����������
�

trainingp"�����������
3__inference_dense_model_7064_layer_call_fn_81179802g;�8
!�
�
x����������
�

trainingp "�����������
3__inference_dense_model_7064_layer_call_fn_81179835g;�8
!�
�
x����������
�

trainingp"�����������
&__inference_signature_wrapper_81179769�<�9
� 
2�/
-
input_1"�
input_1����������"3�0
.
output_1"�
output_1���������