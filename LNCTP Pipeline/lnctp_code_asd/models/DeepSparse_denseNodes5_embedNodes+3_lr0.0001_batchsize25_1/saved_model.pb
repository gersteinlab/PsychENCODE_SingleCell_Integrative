��
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
 �"serve*2.13.02v2.13.0-rc2-7-g1cb1a030a628̈
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
$Adam/v/sparse_model_8/dense_359/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/v/sparse_model_8/dense_359/bias
�
8Adam/v/sparse_model_8/dense_359/bias/Read/ReadVariableOpReadVariableOp$Adam/v/sparse_model_8/dense_359/bias*
_output_shapes
:*
dtype0
�
$Adam/m/sparse_model_8/dense_359/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/m/sparse_model_8/dense_359/bias
�
8Adam/m/sparse_model_8/dense_359/bias/Read/ReadVariableOpReadVariableOp$Adam/m/sparse_model_8/dense_359/bias*
_output_shapes
:*
dtype0
�
&Adam/v/sparse_model_8/dense_359/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/v/sparse_model_8/dense_359/kernel
�
:Adam/v/sparse_model_8/dense_359/kernel/Read/ReadVariableOpReadVariableOp&Adam/v/sparse_model_8/dense_359/kernel*
_output_shapes

:*
dtype0
�
&Adam/m/sparse_model_8/dense_359/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/m/sparse_model_8/dense_359/kernel
�
:Adam/m/sparse_model_8/dense_359/kernel/Read/ReadVariableOpReadVariableOp&Adam/m/sparse_model_8/dense_359/kernel*
_output_shapes

:*
dtype0
�
2Adam/v/sparse_model_8/batch_normalization_242/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/v/sparse_model_8/batch_normalization_242/beta
�
FAdam/v/sparse_model_8/batch_normalization_242/beta/Read/ReadVariableOpReadVariableOp2Adam/v/sparse_model_8/batch_normalization_242/beta*
_output_shapes
:*
dtype0
�
2Adam/m/sparse_model_8/batch_normalization_242/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/m/sparse_model_8/batch_normalization_242/beta
�
FAdam/m/sparse_model_8/batch_normalization_242/beta/Read/ReadVariableOpReadVariableOp2Adam/m/sparse_model_8/batch_normalization_242/beta*
_output_shapes
:*
dtype0
�
3Adam/v/sparse_model_8/batch_normalization_242/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/sparse_model_8/batch_normalization_242/gamma
�
GAdam/v/sparse_model_8/batch_normalization_242/gamma/Read/ReadVariableOpReadVariableOp3Adam/v/sparse_model_8/batch_normalization_242/gamma*
_output_shapes
:*
dtype0
�
3Adam/m/sparse_model_8/batch_normalization_242/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/sparse_model_8/batch_normalization_242/gamma
�
GAdam/m/sparse_model_8/batch_normalization_242/gamma/Read/ReadVariableOpReadVariableOp3Adam/m/sparse_model_8/batch_normalization_242/gamma*
_output_shapes
:*
dtype0
�
$Adam/v/sparse_model_8/dense_358/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/v/sparse_model_8/dense_358/bias
�
8Adam/v/sparse_model_8/dense_358/bias/Read/ReadVariableOpReadVariableOp$Adam/v/sparse_model_8/dense_358/bias*
_output_shapes
:*
dtype0
�
$Adam/m/sparse_model_8/dense_358/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/m/sparse_model_8/dense_358/bias
�
8Adam/m/sparse_model_8/dense_358/bias/Read/ReadVariableOpReadVariableOp$Adam/m/sparse_model_8/dense_358/bias*
_output_shapes
:*
dtype0
�
&Adam/v/sparse_model_8/dense_358/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/v/sparse_model_8/dense_358/kernel
�
:Adam/v/sparse_model_8/dense_358/kernel/Read/ReadVariableOpReadVariableOp&Adam/v/sparse_model_8/dense_358/kernel*
_output_shapes

:*
dtype0
�
&Adam/m/sparse_model_8/dense_358/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/m/sparse_model_8/dense_358/kernel
�
:Adam/m/sparse_model_8/dense_358/kernel/Read/ReadVariableOpReadVariableOp&Adam/m/sparse_model_8/dense_358/kernel*
_output_shapes

:*
dtype0
�
2Adam/v/sparse_model_8/batch_normalization_241/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/v/sparse_model_8/batch_normalization_241/beta
�
FAdam/v/sparse_model_8/batch_normalization_241/beta/Read/ReadVariableOpReadVariableOp2Adam/v/sparse_model_8/batch_normalization_241/beta*
_output_shapes
:*
dtype0
�
2Adam/m/sparse_model_8/batch_normalization_241/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/m/sparse_model_8/batch_normalization_241/beta
�
FAdam/m/sparse_model_8/batch_normalization_241/beta/Read/ReadVariableOpReadVariableOp2Adam/m/sparse_model_8/batch_normalization_241/beta*
_output_shapes
:*
dtype0
�
3Adam/v/sparse_model_8/batch_normalization_241/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/sparse_model_8/batch_normalization_241/gamma
�
GAdam/v/sparse_model_8/batch_normalization_241/gamma/Read/ReadVariableOpReadVariableOp3Adam/v/sparse_model_8/batch_normalization_241/gamma*
_output_shapes
:*
dtype0
�
3Adam/m/sparse_model_8/batch_normalization_241/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/sparse_model_8/batch_normalization_241/gamma
�
GAdam/m/sparse_model_8/batch_normalization_241/gamma/Read/ReadVariableOpReadVariableOp3Adam/m/sparse_model_8/batch_normalization_241/gamma*
_output_shapes
:*
dtype0
�
$Adam/v/sparse_model_8/dense_357/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/v/sparse_model_8/dense_357/bias
�
8Adam/v/sparse_model_8/dense_357/bias/Read/ReadVariableOpReadVariableOp$Adam/v/sparse_model_8/dense_357/bias*
_output_shapes
:*
dtype0
�
$Adam/m/sparse_model_8/dense_357/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/m/sparse_model_8/dense_357/bias
�
8Adam/m/sparse_model_8/dense_357/bias/Read/ReadVariableOpReadVariableOp$Adam/m/sparse_model_8/dense_357/bias*
_output_shapes
:*
dtype0
�
&Adam/v/sparse_model_8/dense_357/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/v/sparse_model_8/dense_357/kernel
�
:Adam/v/sparse_model_8/dense_357/kernel/Read/ReadVariableOpReadVariableOp&Adam/v/sparse_model_8/dense_357/kernel*
_output_shapes

:*
dtype0
�
&Adam/m/sparse_model_8/dense_357/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*7
shared_name(&Adam/m/sparse_model_8/dense_357/kernel
�
:Adam/m/sparse_model_8/dense_357/kernel/Read/ReadVariableOpReadVariableOp&Adam/m/sparse_model_8/dense_357/kernel*
_output_shapes

:*
dtype0
�
2Adam/v/sparse_model_8/batch_normalization_240/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/v/sparse_model_8/batch_normalization_240/beta
�
FAdam/v/sparse_model_8/batch_normalization_240/beta/Read/ReadVariableOpReadVariableOp2Adam/v/sparse_model_8/batch_normalization_240/beta*
_output_shapes
:*
dtype0
�
2Adam/m/sparse_model_8/batch_normalization_240/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42Adam/m/sparse_model_8/batch_normalization_240/beta
�
FAdam/m/sparse_model_8/batch_normalization_240/beta/Read/ReadVariableOpReadVariableOp2Adam/m/sparse_model_8/batch_normalization_240/beta*
_output_shapes
:*
dtype0
�
3Adam/v/sparse_model_8/batch_normalization_240/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/sparse_model_8/batch_normalization_240/gamma
�
GAdam/v/sparse_model_8/batch_normalization_240/gamma/Read/ReadVariableOpReadVariableOp3Adam/v/sparse_model_8/batch_normalization_240/gamma*
_output_shapes
:*
dtype0
�
3Adam/m/sparse_model_8/batch_normalization_240/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/sparse_model_8/batch_normalization_240/gamma
�
GAdam/m/sparse_model_8/batch_normalization_240/gamma/Read/ReadVariableOpReadVariableOp3Adam/m/sparse_model_8/batch_normalization_240/gamma*
_output_shapes
:*
dtype0
�
$Adam/v/sparse_model_8/dense_356/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/v/sparse_model_8/dense_356/bias
�
8Adam/v/sparse_model_8/dense_356/bias/Read/ReadVariableOpReadVariableOp$Adam/v/sparse_model_8/dense_356/bias*
_output_shapes
:*
dtype0
�
$Adam/m/sparse_model_8/dense_356/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*5
shared_name&$Adam/m/sparse_model_8/dense_356/bias
�
8Adam/m/sparse_model_8/dense_356/bias/Read/ReadVariableOpReadVariableOp$Adam/m/sparse_model_8/dense_356/bias*
_output_shapes
:*
dtype0
�
&Adam/v/sparse_model_8/dense_356/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*7
shared_name(&Adam/v/sparse_model_8/dense_356/kernel
�
:Adam/v/sparse_model_8/dense_356/kernel/Read/ReadVariableOpReadVariableOp&Adam/v/sparse_model_8/dense_356/kernel*
_output_shapes
:	�*
dtype0
�
&Adam/m/sparse_model_8/dense_356/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*7
shared_name(&Adam/m/sparse_model_8/dense_356/kernel
�
:Adam/m/sparse_model_8/dense_356/kernel/Read/ReadVariableOpReadVariableOp&Adam/m/sparse_model_8/dense_356/kernel*
_output_shapes
:	�*
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
sparse_model_8/dense_359/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*.
shared_namesparse_model_8/dense_359/bias
�
1sparse_model_8/dense_359/bias/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_359/bias*
_output_shapes
:*
dtype0
�
sparse_model_8/dense_359/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*0
shared_name!sparse_model_8/dense_359/kernel
�
3sparse_model_8/dense_359/kernel/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_359/kernel*
_output_shapes

:*
dtype0
�
6sparse_model_8/batch_normalization_242/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*G
shared_name86sparse_model_8/batch_normalization_242/moving_variance
�
Jsparse_model_8/batch_normalization_242/moving_variance/Read/ReadVariableOpReadVariableOp6sparse_model_8/batch_normalization_242/moving_variance*
_output_shapes
:*
dtype0
�
2sparse_model_8/batch_normalization_242/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42sparse_model_8/batch_normalization_242/moving_mean
�
Fsparse_model_8/batch_normalization_242/moving_mean/Read/ReadVariableOpReadVariableOp2sparse_model_8/batch_normalization_242/moving_mean*
_output_shapes
:*
dtype0
�
+sparse_model_8/batch_normalization_242/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*<
shared_name-+sparse_model_8/batch_normalization_242/beta
�
?sparse_model_8/batch_normalization_242/beta/Read/ReadVariableOpReadVariableOp+sparse_model_8/batch_normalization_242/beta*
_output_shapes
:*
dtype0
�
,sparse_model_8/batch_normalization_242/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,sparse_model_8/batch_normalization_242/gamma
�
@sparse_model_8/batch_normalization_242/gamma/Read/ReadVariableOpReadVariableOp,sparse_model_8/batch_normalization_242/gamma*
_output_shapes
:*
dtype0
�
sparse_model_8/dense_358/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*.
shared_namesparse_model_8/dense_358/bias
�
1sparse_model_8/dense_358/bias/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_358/bias*
_output_shapes
:*
dtype0
�
sparse_model_8/dense_358/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*0
shared_name!sparse_model_8/dense_358/kernel
�
3sparse_model_8/dense_358/kernel/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_358/kernel*
_output_shapes

:*
dtype0
�
6sparse_model_8/batch_normalization_241/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*G
shared_name86sparse_model_8/batch_normalization_241/moving_variance
�
Jsparse_model_8/batch_normalization_241/moving_variance/Read/ReadVariableOpReadVariableOp6sparse_model_8/batch_normalization_241/moving_variance*
_output_shapes
:*
dtype0
�
2sparse_model_8/batch_normalization_241/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42sparse_model_8/batch_normalization_241/moving_mean
�
Fsparse_model_8/batch_normalization_241/moving_mean/Read/ReadVariableOpReadVariableOp2sparse_model_8/batch_normalization_241/moving_mean*
_output_shapes
:*
dtype0
�
+sparse_model_8/batch_normalization_241/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*<
shared_name-+sparse_model_8/batch_normalization_241/beta
�
?sparse_model_8/batch_normalization_241/beta/Read/ReadVariableOpReadVariableOp+sparse_model_8/batch_normalization_241/beta*
_output_shapes
:*
dtype0
�
,sparse_model_8/batch_normalization_241/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,sparse_model_8/batch_normalization_241/gamma
�
@sparse_model_8/batch_normalization_241/gamma/Read/ReadVariableOpReadVariableOp,sparse_model_8/batch_normalization_241/gamma*
_output_shapes
:*
dtype0
�
sparse_model_8/dense_357/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*.
shared_namesparse_model_8/dense_357/bias
�
1sparse_model_8/dense_357/bias/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_357/bias*
_output_shapes
:*
dtype0
�
sparse_model_8/dense_357/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*0
shared_name!sparse_model_8/dense_357/kernel
�
3sparse_model_8/dense_357/kernel/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_357/kernel*
_output_shapes

:*
dtype0
�
6sparse_model_8/batch_normalization_240/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*G
shared_name86sparse_model_8/batch_normalization_240/moving_variance
�
Jsparse_model_8/batch_normalization_240/moving_variance/Read/ReadVariableOpReadVariableOp6sparse_model_8/batch_normalization_240/moving_variance*
_output_shapes
:*
dtype0
�
2sparse_model_8/batch_normalization_240/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*C
shared_name42sparse_model_8/batch_normalization_240/moving_mean
�
Fsparse_model_8/batch_normalization_240/moving_mean/Read/ReadVariableOpReadVariableOp2sparse_model_8/batch_normalization_240/moving_mean*
_output_shapes
:*
dtype0
�
+sparse_model_8/batch_normalization_240/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*<
shared_name-+sparse_model_8/batch_normalization_240/beta
�
?sparse_model_8/batch_normalization_240/beta/Read/ReadVariableOpReadVariableOp+sparse_model_8/batch_normalization_240/beta*
_output_shapes
:*
dtype0
�
,sparse_model_8/batch_normalization_240/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,sparse_model_8/batch_normalization_240/gamma
�
@sparse_model_8/batch_normalization_240/gamma/Read/ReadVariableOpReadVariableOp,sparse_model_8/batch_normalization_240/gamma*
_output_shapes
:*
dtype0
�
sparse_model_8/dense_356/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*.
shared_namesparse_model_8/dense_356/bias
�
1sparse_model_8/dense_356/bias/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_356/bias*
_output_shapes
:*
dtype0
�
sparse_model_8/dense_356/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*0
shared_name!sparse_model_8/dense_356/kernel
�
3sparse_model_8/dense_356/kernel/Read/ReadVariableOpReadVariableOpsparse_model_8/dense_356/kernel*
_output_shapes
:	�*
dtype0
|
serving_default_input_1Placeholder*(
_output_shapes
:����������*
dtype0*
shape:����������
�	
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1sparse_model_8/dense_356/kernelsparse_model_8/dense_356/bias2sparse_model_8/batch_normalization_240/moving_mean6sparse_model_8/batch_normalization_240/moving_variance+sparse_model_8/batch_normalization_240/beta,sparse_model_8/batch_normalization_240/gammasparse_model_8/dense_357/kernelsparse_model_8/dense_357/bias2sparse_model_8/batch_normalization_241/moving_mean6sparse_model_8/batch_normalization_241/moving_variance+sparse_model_8/batch_normalization_241/beta,sparse_model_8/batch_normalization_241/gammasparse_model_8/dense_358/kernelsparse_model_8/dense_358/bias2sparse_model_8/batch_normalization_242/moving_mean6sparse_model_8/batch_normalization_242/moving_variance+sparse_model_8/batch_normalization_242/beta,sparse_model_8/batch_normalization_242/gammasparse_model_8/dense_359/kernelsparse_model_8/dense_359/bias* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *-
f(R&
$__inference_signature_wrapper_448939

NoOpNoOp
�`
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�_
value�_B�_ B�_
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


dense1
	norm2

dense2
	norm3
out
	optimizer

signatures*
�
0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
 15
!16
"17
#18
$19*
j
0
1
2
3
4
5
6
7
8
9
10
 11
#12
$13*
* 
�
%non_trainable_variables

&layers
'metrics
(layer_regularization_losses
)layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*

*trace_0
+trace_1* 

,trace_0
-trace_1* 
* 
�
.	variables
/trainable_variables
0regularization_losses
1	keras_api
2__call__
*3&call_and_return_all_conditional_losses

kernel
bias*
�
4	variables
5trainable_variables
6regularization_losses
7	keras_api
8__call__
*9&call_and_return_all_conditional_losses
:axis
	gamma
beta
moving_mean
moving_variance*
�
;	variables
<trainable_variables
=regularization_losses
>	keras_api
?__call__
*@&call_and_return_all_conditional_losses

kernel
bias*
�
A	variables
Btrainable_variables
Cregularization_losses
D	keras_api
E__call__
*F&call_and_return_all_conditional_losses
Gaxis
	gamma
beta
moving_mean
moving_variance*
�
H	variables
Itrainable_variables
Jregularization_losses
K	keras_api
L__call__
*M&call_and_return_all_conditional_losses

kernel
bias*
�
N	variables
Otrainable_variables
Pregularization_losses
Q	keras_api
R__call__
*S&call_and_return_all_conditional_losses
Taxis
	gamma
 beta
!moving_mean
"moving_variance*
�
U	variables
Vtrainable_variables
Wregularization_losses
X	keras_api
Y__call__
*Z&call_and_return_all_conditional_losses

#kernel
$bias*
�
[
_variables
\_iterations
]_learning_rate
^_index_dict
_
_momentums
`_velocities
a_update_step_xla*

bserving_default* 
_Y
VARIABLE_VALUEsparse_model_8/dense_356/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEsparse_model_8/dense_356/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE,sparse_model_8/batch_normalization_240/gamma&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
ke
VARIABLE_VALUE+sparse_model_8/batch_normalization_240/beta&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE2sparse_model_8/batch_normalization_240/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
vp
VARIABLE_VALUE6sparse_model_8/batch_normalization_240/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEsparse_model_8/dense_357/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEsparse_model_8/dense_357/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE,sparse_model_8/batch_normalization_241/gamma&variables/8/.ATTRIBUTES/VARIABLE_VALUE*
ke
VARIABLE_VALUE+sparse_model_8/batch_normalization_241/beta&variables/9/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE2sparse_model_8/batch_normalization_241/moving_mean'variables/10/.ATTRIBUTES/VARIABLE_VALUE*
wq
VARIABLE_VALUE6sparse_model_8/batch_normalization_241/moving_variance'variables/11/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEsparse_model_8/dense_358/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEsparse_model_8/dense_358/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE*
mg
VARIABLE_VALUE,sparse_model_8/batch_normalization_242/gamma'variables/14/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE+sparse_model_8/batch_normalization_242/beta'variables/15/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE2sparse_model_8/batch_normalization_242/moving_mean'variables/16/.ATTRIBUTES/VARIABLE_VALUE*
wq
VARIABLE_VALUE6sparse_model_8/batch_normalization_242/moving_variance'variables/17/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEsparse_model_8/dense_359/kernel'variables/18/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEsparse_model_8/dense_359/bias'variables/19/.ATTRIBUTES/VARIABLE_VALUE*
.
0
1
2
3
!4
"5*
5
0
	1

2
3
4
5
6*

c0
d1*
* 
* 
* 
* 
* 
* 

0
1*

0
1*
* 
�
enon_trainable_variables

flayers
gmetrics
hlayer_regularization_losses
ilayer_metrics
.	variables
/trainable_variables
0regularization_losses
2__call__
*3&call_and_return_all_conditional_losses
&3"call_and_return_conditional_losses*

jtrace_0* 

ktrace_0* 
 
0
1
2
3*

0
1*
* 
�
lnon_trainable_variables

mlayers
nmetrics
olayer_regularization_losses
player_metrics
4	variables
5trainable_variables
6regularization_losses
8__call__
*9&call_and_return_all_conditional_losses
&9"call_and_return_conditional_losses*

qtrace_0
rtrace_1* 

strace_0
ttrace_1* 
* 

0
1*

0
1*
* 
�
unon_trainable_variables

vlayers
wmetrics
xlayer_regularization_losses
ylayer_metrics
;	variables
<trainable_variables
=regularization_losses
?__call__
*@&call_and_return_all_conditional_losses
&@"call_and_return_conditional_losses*

ztrace_0* 

{trace_0* 
 
0
1
2
3*

0
1*
* 
�
|non_trainable_variables

}layers
~metrics
layer_regularization_losses
�layer_metrics
A	variables
Btrainable_variables
Cregularization_losses
E__call__
*F&call_and_return_all_conditional_losses
&F"call_and_return_conditional_losses*

�trace_0
�trace_1* 

�trace_0
�trace_1* 
* 

0
1*

0
1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
H	variables
Itrainable_variables
Jregularization_losses
L__call__
*M&call_and_return_all_conditional_losses
&M"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
 
0
 1
!2
"3*

0
 1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
N	variables
Otrainable_variables
Pregularization_losses
R__call__
*S&call_and_return_all_conditional_losses
&S"call_and_return_conditional_losses*

�trace_0
�trace_1* 

�trace_0
�trace_1* 
* 

#0
$1*

#0
$1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
U	variables
Vtrainable_variables
Wregularization_losses
Y__call__
*Z&call_and_return_all_conditional_losses
&Z"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
�
\0
�1
�2
�3
�4
�5
�6
�7
�8
�9
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
�20
�21
�22
�23
�24
�25
�26
�27
�28*
SM
VARIABLE_VALUE	iteration0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUElearning_rate3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 
x
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13*
x
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13*
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
0
1*
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
0
1*
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
!0
"1*
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
VARIABLE_VALUE&Adam/m/sparse_model_8/dense_356/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE&Adam/v/sparse_model_8/dense_356/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
oi
VARIABLE_VALUE$Adam/m/sparse_model_8/dense_356/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
oi
VARIABLE_VALUE$Adam/v/sparse_model_8/dense_356/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE3Adam/m/sparse_model_8/batch_normalization_240/gamma1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE3Adam/v/sparse_model_8/batch_normalization_240/gamma1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
}w
VARIABLE_VALUE2Adam/m/sparse_model_8/batch_normalization_240/beta1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
}w
VARIABLE_VALUE2Adam/v/sparse_model_8/batch_normalization_240/beta1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE&Adam/m/sparse_model_8/dense_357/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/v/sparse_model_8/dense_357/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/m/sparse_model_8/dense_357/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/v/sparse_model_8/dense_357/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/m/sparse_model_8/batch_normalization_241/gamma2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/v/sparse_model_8/batch_normalization_241/gamma2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE2Adam/m/sparse_model_8/batch_normalization_241/beta2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE2Adam/v/sparse_model_8/batch_normalization_241/beta2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/m/sparse_model_8/dense_358/kernel2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/v/sparse_model_8/dense_358/kernel2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/m/sparse_model_8/dense_358/bias2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/v/sparse_model_8/dense_358/bias2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/m/sparse_model_8/batch_normalization_242/gamma2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/v/sparse_model_8/batch_normalization_242/gamma2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE2Adam/m/sparse_model_8/batch_normalization_242/beta2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE2Adam/v/sparse_model_8/batch_normalization_242/beta2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/m/sparse_model_8/dense_359/kernel2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE&Adam/v/sparse_model_8/dense_359/kernel2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/m/sparse_model_8/dense_359/bias2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE$Adam/v/sparse_model_8/dense_359/bias2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUE*
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
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenamesparse_model_8/dense_356/kernelsparse_model_8/dense_356/bias,sparse_model_8/batch_normalization_240/gamma+sparse_model_8/batch_normalization_240/beta2sparse_model_8/batch_normalization_240/moving_mean6sparse_model_8/batch_normalization_240/moving_variancesparse_model_8/dense_357/kernelsparse_model_8/dense_357/bias,sparse_model_8/batch_normalization_241/gamma+sparse_model_8/batch_normalization_241/beta2sparse_model_8/batch_normalization_241/moving_mean6sparse_model_8/batch_normalization_241/moving_variancesparse_model_8/dense_358/kernelsparse_model_8/dense_358/bias,sparse_model_8/batch_normalization_242/gamma+sparse_model_8/batch_normalization_242/beta2sparse_model_8/batch_normalization_242/moving_mean6sparse_model_8/batch_normalization_242/moving_variancesparse_model_8/dense_359/kernelsparse_model_8/dense_359/bias	iterationlearning_rate&Adam/m/sparse_model_8/dense_356/kernel&Adam/v/sparse_model_8/dense_356/kernel$Adam/m/sparse_model_8/dense_356/bias$Adam/v/sparse_model_8/dense_356/bias3Adam/m/sparse_model_8/batch_normalization_240/gamma3Adam/v/sparse_model_8/batch_normalization_240/gamma2Adam/m/sparse_model_8/batch_normalization_240/beta2Adam/v/sparse_model_8/batch_normalization_240/beta&Adam/m/sparse_model_8/dense_357/kernel&Adam/v/sparse_model_8/dense_357/kernel$Adam/m/sparse_model_8/dense_357/bias$Adam/v/sparse_model_8/dense_357/bias3Adam/m/sparse_model_8/batch_normalization_241/gamma3Adam/v/sparse_model_8/batch_normalization_241/gamma2Adam/m/sparse_model_8/batch_normalization_241/beta2Adam/v/sparse_model_8/batch_normalization_241/beta&Adam/m/sparse_model_8/dense_358/kernel&Adam/v/sparse_model_8/dense_358/kernel$Adam/m/sparse_model_8/dense_358/bias$Adam/v/sparse_model_8/dense_358/bias3Adam/m/sparse_model_8/batch_normalization_242/gamma3Adam/v/sparse_model_8/batch_normalization_242/gamma2Adam/m/sparse_model_8/batch_normalization_242/beta2Adam/v/sparse_model_8/batch_normalization_242/beta&Adam/m/sparse_model_8/dense_359/kernel&Adam/v/sparse_model_8/dense_359/kernel$Adam/m/sparse_model_8/dense_359/bias$Adam/v/sparse_model_8/dense_359/biastotal_1count_1totalcountConst*C
Tin<
:28*
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
__inference__traced_save_449602
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamesparse_model_8/dense_356/kernelsparse_model_8/dense_356/bias,sparse_model_8/batch_normalization_240/gamma+sparse_model_8/batch_normalization_240/beta2sparse_model_8/batch_normalization_240/moving_mean6sparse_model_8/batch_normalization_240/moving_variancesparse_model_8/dense_357/kernelsparse_model_8/dense_357/bias,sparse_model_8/batch_normalization_241/gamma+sparse_model_8/batch_normalization_241/beta2sparse_model_8/batch_normalization_241/moving_mean6sparse_model_8/batch_normalization_241/moving_variancesparse_model_8/dense_358/kernelsparse_model_8/dense_358/bias,sparse_model_8/batch_normalization_242/gamma+sparse_model_8/batch_normalization_242/beta2sparse_model_8/batch_normalization_242/moving_mean6sparse_model_8/batch_normalization_242/moving_variancesparse_model_8/dense_359/kernelsparse_model_8/dense_359/bias	iterationlearning_rate&Adam/m/sparse_model_8/dense_356/kernel&Adam/v/sparse_model_8/dense_356/kernel$Adam/m/sparse_model_8/dense_356/bias$Adam/v/sparse_model_8/dense_356/bias3Adam/m/sparse_model_8/batch_normalization_240/gamma3Adam/v/sparse_model_8/batch_normalization_240/gamma2Adam/m/sparse_model_8/batch_normalization_240/beta2Adam/v/sparse_model_8/batch_normalization_240/beta&Adam/m/sparse_model_8/dense_357/kernel&Adam/v/sparse_model_8/dense_357/kernel$Adam/m/sparse_model_8/dense_357/bias$Adam/v/sparse_model_8/dense_357/bias3Adam/m/sparse_model_8/batch_normalization_241/gamma3Adam/v/sparse_model_8/batch_normalization_241/gamma2Adam/m/sparse_model_8/batch_normalization_241/beta2Adam/v/sparse_model_8/batch_normalization_241/beta&Adam/m/sparse_model_8/dense_358/kernel&Adam/v/sparse_model_8/dense_358/kernel$Adam/m/sparse_model_8/dense_358/bias$Adam/v/sparse_model_8/dense_358/bias3Adam/m/sparse_model_8/batch_normalization_242/gamma3Adam/v/sparse_model_8/batch_normalization_242/gamma2Adam/m/sparse_model_8/batch_normalization_242/beta2Adam/v/sparse_model_8/batch_normalization_242/beta&Adam/m/sparse_model_8/dense_359/kernel&Adam/v/sparse_model_8/dense_359/kernel$Adam/m/sparse_model_8/dense_359/bias$Adam/v/sparse_model_8/dense_359/biastotal_1count_1totalcount*B
Tin;
927*
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
"__inference__traced_restore_449773��
�	
�
E__inference_dense_358_layer_call_and_return_conditional_losses_448682

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
E__inference_dense_356_layer_call_and_return_conditional_losses_448634

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
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
8__inference_batch_normalization_242_layer_call_fn_449182

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
GPU 2J 8� *\
fWRU
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_448596o
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
 
_user_specified_name449178:&"
 
_user_specified_name449176:&"
 
_user_specified_name449174:&"
 
_user_specified_name449172:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
8__inference_batch_normalization_241_layer_call_fn_449070

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
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_448496o
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
 
_user_specified_name449066:&"
 
_user_specified_name449064:&"
 
_user_specified_name449062:&"
 
_user_specified_name449060:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_448516

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
��
�8
__inference__traced_save_449602
file_prefixI
6read_disablecopyonread_sparse_model_8_dense_356_kernel:	�D
6read_1_disablecopyonread_sparse_model_8_dense_356_bias:S
Eread_2_disablecopyonread_sparse_model_8_batch_normalization_240_gamma:R
Dread_3_disablecopyonread_sparse_model_8_batch_normalization_240_beta:Y
Kread_4_disablecopyonread_sparse_model_8_batch_normalization_240_moving_mean:]
Oread_5_disablecopyonread_sparse_model_8_batch_normalization_240_moving_variance:J
8read_6_disablecopyonread_sparse_model_8_dense_357_kernel:D
6read_7_disablecopyonread_sparse_model_8_dense_357_bias:S
Eread_8_disablecopyonread_sparse_model_8_batch_normalization_241_gamma:R
Dread_9_disablecopyonread_sparse_model_8_batch_normalization_241_beta:Z
Lread_10_disablecopyonread_sparse_model_8_batch_normalization_241_moving_mean:^
Pread_11_disablecopyonread_sparse_model_8_batch_normalization_241_moving_variance:K
9read_12_disablecopyonread_sparse_model_8_dense_358_kernel:E
7read_13_disablecopyonread_sparse_model_8_dense_358_bias:T
Fread_14_disablecopyonread_sparse_model_8_batch_normalization_242_gamma:S
Eread_15_disablecopyonread_sparse_model_8_batch_normalization_242_beta:Z
Lread_16_disablecopyonread_sparse_model_8_batch_normalization_242_moving_mean:^
Pread_17_disablecopyonread_sparse_model_8_batch_normalization_242_moving_variance:K
9read_18_disablecopyonread_sparse_model_8_dense_359_kernel:E
7read_19_disablecopyonread_sparse_model_8_dense_359_bias:-
#read_20_disablecopyonread_iteration:	 1
'read_21_disablecopyonread_learning_rate: S
@read_22_disablecopyonread_adam_m_sparse_model_8_dense_356_kernel:	�S
@read_23_disablecopyonread_adam_v_sparse_model_8_dense_356_kernel:	�L
>read_24_disablecopyonread_adam_m_sparse_model_8_dense_356_bias:L
>read_25_disablecopyonread_adam_v_sparse_model_8_dense_356_bias:[
Mread_26_disablecopyonread_adam_m_sparse_model_8_batch_normalization_240_gamma:[
Mread_27_disablecopyonread_adam_v_sparse_model_8_batch_normalization_240_gamma:Z
Lread_28_disablecopyonread_adam_m_sparse_model_8_batch_normalization_240_beta:Z
Lread_29_disablecopyonread_adam_v_sparse_model_8_batch_normalization_240_beta:R
@read_30_disablecopyonread_adam_m_sparse_model_8_dense_357_kernel:R
@read_31_disablecopyonread_adam_v_sparse_model_8_dense_357_kernel:L
>read_32_disablecopyonread_adam_m_sparse_model_8_dense_357_bias:L
>read_33_disablecopyonread_adam_v_sparse_model_8_dense_357_bias:[
Mread_34_disablecopyonread_adam_m_sparse_model_8_batch_normalization_241_gamma:[
Mread_35_disablecopyonread_adam_v_sparse_model_8_batch_normalization_241_gamma:Z
Lread_36_disablecopyonread_adam_m_sparse_model_8_batch_normalization_241_beta:Z
Lread_37_disablecopyonread_adam_v_sparse_model_8_batch_normalization_241_beta:R
@read_38_disablecopyonread_adam_m_sparse_model_8_dense_358_kernel:R
@read_39_disablecopyonread_adam_v_sparse_model_8_dense_358_kernel:L
>read_40_disablecopyonread_adam_m_sparse_model_8_dense_358_bias:L
>read_41_disablecopyonread_adam_v_sparse_model_8_dense_358_bias:[
Mread_42_disablecopyonread_adam_m_sparse_model_8_batch_normalization_242_gamma:[
Mread_43_disablecopyonread_adam_v_sparse_model_8_batch_normalization_242_gamma:Z
Lread_44_disablecopyonread_adam_m_sparse_model_8_batch_normalization_242_beta:Z
Lread_45_disablecopyonread_adam_v_sparse_model_8_batch_normalization_242_beta:R
@read_46_disablecopyonread_adam_m_sparse_model_8_dense_359_kernel:R
@read_47_disablecopyonread_adam_v_sparse_model_8_dense_359_kernel:L
>read_48_disablecopyonread_adam_m_sparse_model_8_dense_359_bias:L
>read_49_disablecopyonread_adam_v_sparse_model_8_dense_359_bias:+
!read_50_disablecopyonread_total_1: +
!read_51_disablecopyonread_count_1: )
read_52_disablecopyonread_total: )
read_53_disablecopyonread_count: 
savev2_const
identity_109��MergeV2Checkpoints�Read/DisableCopyOnRead�Read/ReadVariableOp�Read_1/DisableCopyOnRead�Read_1/ReadVariableOp�Read_10/DisableCopyOnRead�Read_10/ReadVariableOp�Read_11/DisableCopyOnRead�Read_11/ReadVariableOp�Read_12/DisableCopyOnRead�Read_12/ReadVariableOp�Read_13/DisableCopyOnRead�Read_13/ReadVariableOp�Read_14/DisableCopyOnRead�Read_14/ReadVariableOp�Read_15/DisableCopyOnRead�Read_15/ReadVariableOp�Read_16/DisableCopyOnRead�Read_16/ReadVariableOp�Read_17/DisableCopyOnRead�Read_17/ReadVariableOp�Read_18/DisableCopyOnRead�Read_18/ReadVariableOp�Read_19/DisableCopyOnRead�Read_19/ReadVariableOp�Read_2/DisableCopyOnRead�Read_2/ReadVariableOp�Read_20/DisableCopyOnRead�Read_20/ReadVariableOp�Read_21/DisableCopyOnRead�Read_21/ReadVariableOp�Read_22/DisableCopyOnRead�Read_22/ReadVariableOp�Read_23/DisableCopyOnRead�Read_23/ReadVariableOp�Read_24/DisableCopyOnRead�Read_24/ReadVariableOp�Read_25/DisableCopyOnRead�Read_25/ReadVariableOp�Read_26/DisableCopyOnRead�Read_26/ReadVariableOp�Read_27/DisableCopyOnRead�Read_27/ReadVariableOp�Read_28/DisableCopyOnRead�Read_28/ReadVariableOp�Read_29/DisableCopyOnRead�Read_29/ReadVariableOp�Read_3/DisableCopyOnRead�Read_3/ReadVariableOp�Read_30/DisableCopyOnRead�Read_30/ReadVariableOp�Read_31/DisableCopyOnRead�Read_31/ReadVariableOp�Read_32/DisableCopyOnRead�Read_32/ReadVariableOp�Read_33/DisableCopyOnRead�Read_33/ReadVariableOp�Read_34/DisableCopyOnRead�Read_34/ReadVariableOp�Read_35/DisableCopyOnRead�Read_35/ReadVariableOp�Read_36/DisableCopyOnRead�Read_36/ReadVariableOp�Read_37/DisableCopyOnRead�Read_37/ReadVariableOp�Read_38/DisableCopyOnRead�Read_38/ReadVariableOp�Read_39/DisableCopyOnRead�Read_39/ReadVariableOp�Read_4/DisableCopyOnRead�Read_4/ReadVariableOp�Read_40/DisableCopyOnRead�Read_40/ReadVariableOp�Read_41/DisableCopyOnRead�Read_41/ReadVariableOp�Read_42/DisableCopyOnRead�Read_42/ReadVariableOp�Read_43/DisableCopyOnRead�Read_43/ReadVariableOp�Read_44/DisableCopyOnRead�Read_44/ReadVariableOp�Read_45/DisableCopyOnRead�Read_45/ReadVariableOp�Read_46/DisableCopyOnRead�Read_46/ReadVariableOp�Read_47/DisableCopyOnRead�Read_47/ReadVariableOp�Read_48/DisableCopyOnRead�Read_48/ReadVariableOp�Read_49/DisableCopyOnRead�Read_49/ReadVariableOp�Read_5/DisableCopyOnRead�Read_5/ReadVariableOp�Read_50/DisableCopyOnRead�Read_50/ReadVariableOp�Read_51/DisableCopyOnRead�Read_51/ReadVariableOp�Read_52/DisableCopyOnRead�Read_52/ReadVariableOp�Read_53/DisableCopyOnRead�Read_53/ReadVariableOp�Read_6/DisableCopyOnRead�Read_6/ReadVariableOp�Read_7/DisableCopyOnRead�Read_7/ReadVariableOp�Read_8/DisableCopyOnRead�Read_8/ReadVariableOp�Read_9/DisableCopyOnRead�Read_9/ReadVariableOpw
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
Read/DisableCopyOnReadDisableCopyOnRead6read_disablecopyonread_sparse_model_8_dense_356_kernel"/device:CPU:0*
_output_shapes
 �
Read/ReadVariableOpReadVariableOp6read_disablecopyonread_sparse_model_8_dense_356_kernel^Read/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0j
IdentityIdentityRead/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�b

Identity_1IdentityIdentity:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_1/DisableCopyOnReadDisableCopyOnRead6read_1_disablecopyonread_sparse_model_8_dense_356_bias"/device:CPU:0*
_output_shapes
 �
Read_1/ReadVariableOpReadVariableOp6read_1_disablecopyonread_sparse_model_8_dense_356_bias^Read_1/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_2IdentityRead_1/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_3IdentityIdentity_2:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_2/DisableCopyOnReadDisableCopyOnReadEread_2_disablecopyonread_sparse_model_8_batch_normalization_240_gamma"/device:CPU:0*
_output_shapes
 �
Read_2/ReadVariableOpReadVariableOpEread_2_disablecopyonread_sparse_model_8_batch_normalization_240_gamma^Read_2/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_4IdentityRead_2/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_5IdentityIdentity_4:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_3/DisableCopyOnReadDisableCopyOnReadDread_3_disablecopyonread_sparse_model_8_batch_normalization_240_beta"/device:CPU:0*
_output_shapes
 �
Read_3/ReadVariableOpReadVariableOpDread_3_disablecopyonread_sparse_model_8_batch_normalization_240_beta^Read_3/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_6IdentityRead_3/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_7IdentityIdentity_6:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_4/DisableCopyOnReadDisableCopyOnReadKread_4_disablecopyonread_sparse_model_8_batch_normalization_240_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_4/ReadVariableOpReadVariableOpKread_4_disablecopyonread_sparse_model_8_batch_normalization_240_moving_mean^Read_4/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0i

Identity_8IdentityRead_4/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:_

Identity_9IdentityIdentity_8:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_5/DisableCopyOnReadDisableCopyOnReadOread_5_disablecopyonread_sparse_model_8_batch_normalization_240_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_5/ReadVariableOpReadVariableOpOread_5_disablecopyonread_sparse_model_8_batch_normalization_240_moving_variance^Read_5/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0j
Identity_10IdentityRead_5/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_11IdentityIdentity_10:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_6/DisableCopyOnReadDisableCopyOnRead8read_6_disablecopyonread_sparse_model_8_dense_357_kernel"/device:CPU:0*
_output_shapes
 �
Read_6/ReadVariableOpReadVariableOp8read_6_disablecopyonread_sparse_model_8_dense_357_kernel^Read_6/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0n
Identity_12IdentityRead_6/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_13IdentityIdentity_12:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_7/DisableCopyOnReadDisableCopyOnRead6read_7_disablecopyonread_sparse_model_8_dense_357_bias"/device:CPU:0*
_output_shapes
 �
Read_7/ReadVariableOpReadVariableOp6read_7_disablecopyonread_sparse_model_8_dense_357_bias^Read_7/DisableCopyOnRead"/device:CPU:0*
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
Read_8/DisableCopyOnReadDisableCopyOnReadEread_8_disablecopyonread_sparse_model_8_batch_normalization_241_gamma"/device:CPU:0*
_output_shapes
 �
Read_8/ReadVariableOpReadVariableOpEread_8_disablecopyonread_sparse_model_8_batch_normalization_241_gamma^Read_8/DisableCopyOnRead"/device:CPU:0*
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
Read_9/DisableCopyOnReadDisableCopyOnReadDread_9_disablecopyonread_sparse_model_8_batch_normalization_241_beta"/device:CPU:0*
_output_shapes
 �
Read_9/ReadVariableOpReadVariableOpDread_9_disablecopyonread_sparse_model_8_batch_normalization_241_beta^Read_9/DisableCopyOnRead"/device:CPU:0*
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
Read_10/DisableCopyOnReadDisableCopyOnReadLread_10_disablecopyonread_sparse_model_8_batch_normalization_241_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_10/ReadVariableOpReadVariableOpLread_10_disablecopyonread_sparse_model_8_batch_normalization_241_moving_mean^Read_10/DisableCopyOnRead"/device:CPU:0*
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
Read_11/DisableCopyOnReadDisableCopyOnReadPread_11_disablecopyonread_sparse_model_8_batch_normalization_241_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_11/ReadVariableOpReadVariableOpPread_11_disablecopyonread_sparse_model_8_batch_normalization_241_moving_variance^Read_11/DisableCopyOnRead"/device:CPU:0*
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
Read_12/DisableCopyOnReadDisableCopyOnRead9read_12_disablecopyonread_sparse_model_8_dense_358_kernel"/device:CPU:0*
_output_shapes
 �
Read_12/ReadVariableOpReadVariableOp9read_12_disablecopyonread_sparse_model_8_dense_358_kernel^Read_12/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_24IdentityRead_12/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_25IdentityIdentity_24:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_13/DisableCopyOnReadDisableCopyOnRead7read_13_disablecopyonread_sparse_model_8_dense_358_bias"/device:CPU:0*
_output_shapes
 �
Read_13/ReadVariableOpReadVariableOp7read_13_disablecopyonread_sparse_model_8_dense_358_bias^Read_13/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_26IdentityRead_13/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_27IdentityIdentity_26:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_14/DisableCopyOnReadDisableCopyOnReadFread_14_disablecopyonread_sparse_model_8_batch_normalization_242_gamma"/device:CPU:0*
_output_shapes
 �
Read_14/ReadVariableOpReadVariableOpFread_14_disablecopyonread_sparse_model_8_batch_normalization_242_gamma^Read_14/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_28IdentityRead_14/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_29IdentityIdentity_28:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_15/DisableCopyOnReadDisableCopyOnReadEread_15_disablecopyonread_sparse_model_8_batch_normalization_242_beta"/device:CPU:0*
_output_shapes
 �
Read_15/ReadVariableOpReadVariableOpEread_15_disablecopyonread_sparse_model_8_batch_normalization_242_beta^Read_15/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_30IdentityRead_15/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_31IdentityIdentity_30:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_16/DisableCopyOnReadDisableCopyOnReadLread_16_disablecopyonread_sparse_model_8_batch_normalization_242_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_16/ReadVariableOpReadVariableOpLread_16_disablecopyonread_sparse_model_8_batch_normalization_242_moving_mean^Read_16/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_32IdentityRead_16/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_33IdentityIdentity_32:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_17/DisableCopyOnReadDisableCopyOnReadPread_17_disablecopyonread_sparse_model_8_batch_normalization_242_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_17/ReadVariableOpReadVariableOpPread_17_disablecopyonread_sparse_model_8_batch_normalization_242_moving_variance^Read_17/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_34IdentityRead_17/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_35IdentityIdentity_34:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_18/DisableCopyOnReadDisableCopyOnRead9read_18_disablecopyonread_sparse_model_8_dense_359_kernel"/device:CPU:0*
_output_shapes
 �
Read_18/ReadVariableOpReadVariableOp9read_18_disablecopyonread_sparse_model_8_dense_359_kernel^Read_18/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_36IdentityRead_18/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_37IdentityIdentity_36:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_19/DisableCopyOnReadDisableCopyOnRead7read_19_disablecopyonread_sparse_model_8_dense_359_bias"/device:CPU:0*
_output_shapes
 �
Read_19/ReadVariableOpReadVariableOp7read_19_disablecopyonread_sparse_model_8_dense_359_bias^Read_19/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_38IdentityRead_19/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_39IdentityIdentity_38:output:0"/device:CPU:0*
T0*
_output_shapes
:x
Read_20/DisableCopyOnReadDisableCopyOnRead#read_20_disablecopyonread_iteration"/device:CPU:0*
_output_shapes
 �
Read_20/ReadVariableOpReadVariableOp#read_20_disablecopyonread_iteration^Read_20/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0	g
Identity_40IdentityRead_20/ReadVariableOp:value:0"/device:CPU:0*
T0	*
_output_shapes
: ]
Identity_41IdentityIdentity_40:output:0"/device:CPU:0*
T0	*
_output_shapes
: |
Read_21/DisableCopyOnReadDisableCopyOnRead'read_21_disablecopyonread_learning_rate"/device:CPU:0*
_output_shapes
 �
Read_21/ReadVariableOpReadVariableOp'read_21_disablecopyonread_learning_rate^Read_21/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_42IdentityRead_21/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_43IdentityIdentity_42:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_22/DisableCopyOnReadDisableCopyOnRead@read_22_disablecopyonread_adam_m_sparse_model_8_dense_356_kernel"/device:CPU:0*
_output_shapes
 �
Read_22/ReadVariableOpReadVariableOp@read_22_disablecopyonread_adam_m_sparse_model_8_dense_356_kernel^Read_22/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0p
Identity_44IdentityRead_22/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�f
Identity_45IdentityIdentity_44:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_23/DisableCopyOnReadDisableCopyOnRead@read_23_disablecopyonread_adam_v_sparse_model_8_dense_356_kernel"/device:CPU:0*
_output_shapes
 �
Read_23/ReadVariableOpReadVariableOp@read_23_disablecopyonread_adam_v_sparse_model_8_dense_356_kernel^Read_23/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�*
dtype0p
Identity_46IdentityRead_23/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�f
Identity_47IdentityIdentity_46:output:0"/device:CPU:0*
T0*
_output_shapes
:	��
Read_24/DisableCopyOnReadDisableCopyOnRead>read_24_disablecopyonread_adam_m_sparse_model_8_dense_356_bias"/device:CPU:0*
_output_shapes
 �
Read_24/ReadVariableOpReadVariableOp>read_24_disablecopyonread_adam_m_sparse_model_8_dense_356_bias^Read_24/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_48IdentityRead_24/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_49IdentityIdentity_48:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_25/DisableCopyOnReadDisableCopyOnRead>read_25_disablecopyonread_adam_v_sparse_model_8_dense_356_bias"/device:CPU:0*
_output_shapes
 �
Read_25/ReadVariableOpReadVariableOp>read_25_disablecopyonread_adam_v_sparse_model_8_dense_356_bias^Read_25/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_50IdentityRead_25/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_51IdentityIdentity_50:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_26/DisableCopyOnReadDisableCopyOnReadMread_26_disablecopyonread_adam_m_sparse_model_8_batch_normalization_240_gamma"/device:CPU:0*
_output_shapes
 �
Read_26/ReadVariableOpReadVariableOpMread_26_disablecopyonread_adam_m_sparse_model_8_batch_normalization_240_gamma^Read_26/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_52IdentityRead_26/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_53IdentityIdentity_52:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_27/DisableCopyOnReadDisableCopyOnReadMread_27_disablecopyonread_adam_v_sparse_model_8_batch_normalization_240_gamma"/device:CPU:0*
_output_shapes
 �
Read_27/ReadVariableOpReadVariableOpMread_27_disablecopyonread_adam_v_sparse_model_8_batch_normalization_240_gamma^Read_27/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_54IdentityRead_27/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_55IdentityIdentity_54:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_28/DisableCopyOnReadDisableCopyOnReadLread_28_disablecopyonread_adam_m_sparse_model_8_batch_normalization_240_beta"/device:CPU:0*
_output_shapes
 �
Read_28/ReadVariableOpReadVariableOpLread_28_disablecopyonread_adam_m_sparse_model_8_batch_normalization_240_beta^Read_28/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_56IdentityRead_28/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_57IdentityIdentity_56:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_29/DisableCopyOnReadDisableCopyOnReadLread_29_disablecopyonread_adam_v_sparse_model_8_batch_normalization_240_beta"/device:CPU:0*
_output_shapes
 �
Read_29/ReadVariableOpReadVariableOpLread_29_disablecopyonread_adam_v_sparse_model_8_batch_normalization_240_beta^Read_29/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_58IdentityRead_29/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_59IdentityIdentity_58:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_30/DisableCopyOnReadDisableCopyOnRead@read_30_disablecopyonread_adam_m_sparse_model_8_dense_357_kernel"/device:CPU:0*
_output_shapes
 �
Read_30/ReadVariableOpReadVariableOp@read_30_disablecopyonread_adam_m_sparse_model_8_dense_357_kernel^Read_30/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_60IdentityRead_30/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_61IdentityIdentity_60:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_31/DisableCopyOnReadDisableCopyOnRead@read_31_disablecopyonread_adam_v_sparse_model_8_dense_357_kernel"/device:CPU:0*
_output_shapes
 �
Read_31/ReadVariableOpReadVariableOp@read_31_disablecopyonread_adam_v_sparse_model_8_dense_357_kernel^Read_31/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_62IdentityRead_31/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_63IdentityIdentity_62:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_32/DisableCopyOnReadDisableCopyOnRead>read_32_disablecopyonread_adam_m_sparse_model_8_dense_357_bias"/device:CPU:0*
_output_shapes
 �
Read_32/ReadVariableOpReadVariableOp>read_32_disablecopyonread_adam_m_sparse_model_8_dense_357_bias^Read_32/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_64IdentityRead_32/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_65IdentityIdentity_64:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_33/DisableCopyOnReadDisableCopyOnRead>read_33_disablecopyonread_adam_v_sparse_model_8_dense_357_bias"/device:CPU:0*
_output_shapes
 �
Read_33/ReadVariableOpReadVariableOp>read_33_disablecopyonread_adam_v_sparse_model_8_dense_357_bias^Read_33/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_66IdentityRead_33/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_67IdentityIdentity_66:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_34/DisableCopyOnReadDisableCopyOnReadMread_34_disablecopyonread_adam_m_sparse_model_8_batch_normalization_241_gamma"/device:CPU:0*
_output_shapes
 �
Read_34/ReadVariableOpReadVariableOpMread_34_disablecopyonread_adam_m_sparse_model_8_batch_normalization_241_gamma^Read_34/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_68IdentityRead_34/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_69IdentityIdentity_68:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_35/DisableCopyOnReadDisableCopyOnReadMread_35_disablecopyonread_adam_v_sparse_model_8_batch_normalization_241_gamma"/device:CPU:0*
_output_shapes
 �
Read_35/ReadVariableOpReadVariableOpMread_35_disablecopyonread_adam_v_sparse_model_8_batch_normalization_241_gamma^Read_35/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_70IdentityRead_35/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_71IdentityIdentity_70:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_36/DisableCopyOnReadDisableCopyOnReadLread_36_disablecopyonread_adam_m_sparse_model_8_batch_normalization_241_beta"/device:CPU:0*
_output_shapes
 �
Read_36/ReadVariableOpReadVariableOpLread_36_disablecopyonread_adam_m_sparse_model_8_batch_normalization_241_beta^Read_36/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_72IdentityRead_36/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_73IdentityIdentity_72:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_37/DisableCopyOnReadDisableCopyOnReadLread_37_disablecopyonread_adam_v_sparse_model_8_batch_normalization_241_beta"/device:CPU:0*
_output_shapes
 �
Read_37/ReadVariableOpReadVariableOpLread_37_disablecopyonread_adam_v_sparse_model_8_batch_normalization_241_beta^Read_37/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_74IdentityRead_37/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_75IdentityIdentity_74:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_38/DisableCopyOnReadDisableCopyOnRead@read_38_disablecopyonread_adam_m_sparse_model_8_dense_358_kernel"/device:CPU:0*
_output_shapes
 �
Read_38/ReadVariableOpReadVariableOp@read_38_disablecopyonread_adam_m_sparse_model_8_dense_358_kernel^Read_38/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_76IdentityRead_38/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_77IdentityIdentity_76:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_39/DisableCopyOnReadDisableCopyOnRead@read_39_disablecopyonread_adam_v_sparse_model_8_dense_358_kernel"/device:CPU:0*
_output_shapes
 �
Read_39/ReadVariableOpReadVariableOp@read_39_disablecopyonread_adam_v_sparse_model_8_dense_358_kernel^Read_39/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_78IdentityRead_39/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_79IdentityIdentity_78:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_40/DisableCopyOnReadDisableCopyOnRead>read_40_disablecopyonread_adam_m_sparse_model_8_dense_358_bias"/device:CPU:0*
_output_shapes
 �
Read_40/ReadVariableOpReadVariableOp>read_40_disablecopyonread_adam_m_sparse_model_8_dense_358_bias^Read_40/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_80IdentityRead_40/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_81IdentityIdentity_80:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_41/DisableCopyOnReadDisableCopyOnRead>read_41_disablecopyonread_adam_v_sparse_model_8_dense_358_bias"/device:CPU:0*
_output_shapes
 �
Read_41/ReadVariableOpReadVariableOp>read_41_disablecopyonread_adam_v_sparse_model_8_dense_358_bias^Read_41/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_82IdentityRead_41/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_83IdentityIdentity_82:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_42/DisableCopyOnReadDisableCopyOnReadMread_42_disablecopyonread_adam_m_sparse_model_8_batch_normalization_242_gamma"/device:CPU:0*
_output_shapes
 �
Read_42/ReadVariableOpReadVariableOpMread_42_disablecopyonread_adam_m_sparse_model_8_batch_normalization_242_gamma^Read_42/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_84IdentityRead_42/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_85IdentityIdentity_84:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_43/DisableCopyOnReadDisableCopyOnReadMread_43_disablecopyonread_adam_v_sparse_model_8_batch_normalization_242_gamma"/device:CPU:0*
_output_shapes
 �
Read_43/ReadVariableOpReadVariableOpMread_43_disablecopyonread_adam_v_sparse_model_8_batch_normalization_242_gamma^Read_43/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_86IdentityRead_43/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_87IdentityIdentity_86:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_44/DisableCopyOnReadDisableCopyOnReadLread_44_disablecopyonread_adam_m_sparse_model_8_batch_normalization_242_beta"/device:CPU:0*
_output_shapes
 �
Read_44/ReadVariableOpReadVariableOpLread_44_disablecopyonread_adam_m_sparse_model_8_batch_normalization_242_beta^Read_44/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_88IdentityRead_44/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_89IdentityIdentity_88:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_45/DisableCopyOnReadDisableCopyOnReadLread_45_disablecopyonread_adam_v_sparse_model_8_batch_normalization_242_beta"/device:CPU:0*
_output_shapes
 �
Read_45/ReadVariableOpReadVariableOpLread_45_disablecopyonread_adam_v_sparse_model_8_batch_normalization_242_beta^Read_45/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_90IdentityRead_45/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_91IdentityIdentity_90:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_46/DisableCopyOnReadDisableCopyOnRead@read_46_disablecopyonread_adam_m_sparse_model_8_dense_359_kernel"/device:CPU:0*
_output_shapes
 �
Read_46/ReadVariableOpReadVariableOp@read_46_disablecopyonread_adam_m_sparse_model_8_dense_359_kernel^Read_46/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_92IdentityRead_46/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_93IdentityIdentity_92:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_47/DisableCopyOnReadDisableCopyOnRead@read_47_disablecopyonread_adam_v_sparse_model_8_dense_359_kernel"/device:CPU:0*
_output_shapes
 �
Read_47/ReadVariableOpReadVariableOp@read_47_disablecopyonread_adam_v_sparse_model_8_dense_359_kernel^Read_47/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_94IdentityRead_47/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_95IdentityIdentity_94:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_48/DisableCopyOnReadDisableCopyOnRead>read_48_disablecopyonread_adam_m_sparse_model_8_dense_359_bias"/device:CPU:0*
_output_shapes
 �
Read_48/ReadVariableOpReadVariableOp>read_48_disablecopyonread_adam_m_sparse_model_8_dense_359_bias^Read_48/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_96IdentityRead_48/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_97IdentityIdentity_96:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_49/DisableCopyOnReadDisableCopyOnRead>read_49_disablecopyonread_adam_v_sparse_model_8_dense_359_bias"/device:CPU:0*
_output_shapes
 �
Read_49/ReadVariableOpReadVariableOp>read_49_disablecopyonread_adam_v_sparse_model_8_dense_359_bias^Read_49/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_98IdentityRead_49/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_99IdentityIdentity_98:output:0"/device:CPU:0*
T0*
_output_shapes
:v
Read_50/DisableCopyOnReadDisableCopyOnRead!read_50_disablecopyonread_total_1"/device:CPU:0*
_output_shapes
 �
Read_50/ReadVariableOpReadVariableOp!read_50_disablecopyonread_total_1^Read_50/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0h
Identity_100IdentityRead_50/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: _
Identity_101IdentityIdentity_100:output:0"/device:CPU:0*
T0*
_output_shapes
: v
Read_51/DisableCopyOnReadDisableCopyOnRead!read_51_disablecopyonread_count_1"/device:CPU:0*
_output_shapes
 �
Read_51/ReadVariableOpReadVariableOp!read_51_disablecopyonread_count_1^Read_51/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0h
Identity_102IdentityRead_51/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: _
Identity_103IdentityIdentity_102:output:0"/device:CPU:0*
T0*
_output_shapes
: t
Read_52/DisableCopyOnReadDisableCopyOnReadread_52_disablecopyonread_total"/device:CPU:0*
_output_shapes
 �
Read_52/ReadVariableOpReadVariableOpread_52_disablecopyonread_total^Read_52/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0h
Identity_104IdentityRead_52/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: _
Identity_105IdentityIdentity_104:output:0"/device:CPU:0*
T0*
_output_shapes
: t
Read_53/DisableCopyOnReadDisableCopyOnReadread_53_disablecopyonread_count"/device:CPU:0*
_output_shapes
 �
Read_53/ReadVariableOpReadVariableOpread_53_disablecopyonread_count^Read_53/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0h
Identity_106IdentityRead_53/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: _
Identity_107IdentityIdentity_106:output:0"/device:CPU:0*
T0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:7*
dtype0*�
value�B�7B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB'variables/18/.ATTRIBUTES/VARIABLE_VALUEB'variables/19/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:7*
dtype0*�
valuexBv7B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0Identity_1:output:0Identity_3:output:0Identity_5:output:0Identity_7:output:0Identity_9:output:0Identity_11:output:0Identity_13:output:0Identity_15:output:0Identity_17:output:0Identity_19:output:0Identity_21:output:0Identity_23:output:0Identity_25:output:0Identity_27:output:0Identity_29:output:0Identity_31:output:0Identity_33:output:0Identity_35:output:0Identity_37:output:0Identity_39:output:0Identity_41:output:0Identity_43:output:0Identity_45:output:0Identity_47:output:0Identity_49:output:0Identity_51:output:0Identity_53:output:0Identity_55:output:0Identity_57:output:0Identity_59:output:0Identity_61:output:0Identity_63:output:0Identity_65:output:0Identity_67:output:0Identity_69:output:0Identity_71:output:0Identity_73:output:0Identity_75:output:0Identity_77:output:0Identity_79:output:0Identity_81:output:0Identity_83:output:0Identity_85:output:0Identity_87:output:0Identity_89:output:0Identity_91:output:0Identity_93:output:0Identity_95:output:0Identity_97:output:0Identity_99:output:0Identity_101:output:0Identity_103:output:0Identity_105:output:0Identity_107:output:0savev2_const"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *E
dtypes;
927	�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 j
Identity_108Identityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: W
Identity_109IdentityIdentity_108:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^MergeV2Checkpoints^Read/DisableCopyOnRead^Read/ReadVariableOp^Read_1/DisableCopyOnRead^Read_1/ReadVariableOp^Read_10/DisableCopyOnRead^Read_10/ReadVariableOp^Read_11/DisableCopyOnRead^Read_11/ReadVariableOp^Read_12/DisableCopyOnRead^Read_12/ReadVariableOp^Read_13/DisableCopyOnRead^Read_13/ReadVariableOp^Read_14/DisableCopyOnRead^Read_14/ReadVariableOp^Read_15/DisableCopyOnRead^Read_15/ReadVariableOp^Read_16/DisableCopyOnRead^Read_16/ReadVariableOp^Read_17/DisableCopyOnRead^Read_17/ReadVariableOp^Read_18/DisableCopyOnRead^Read_18/ReadVariableOp^Read_19/DisableCopyOnRead^Read_19/ReadVariableOp^Read_2/DisableCopyOnRead^Read_2/ReadVariableOp^Read_20/DisableCopyOnRead^Read_20/ReadVariableOp^Read_21/DisableCopyOnRead^Read_21/ReadVariableOp^Read_22/DisableCopyOnRead^Read_22/ReadVariableOp^Read_23/DisableCopyOnRead^Read_23/ReadVariableOp^Read_24/DisableCopyOnRead^Read_24/ReadVariableOp^Read_25/DisableCopyOnRead^Read_25/ReadVariableOp^Read_26/DisableCopyOnRead^Read_26/ReadVariableOp^Read_27/DisableCopyOnRead^Read_27/ReadVariableOp^Read_28/DisableCopyOnRead^Read_28/ReadVariableOp^Read_29/DisableCopyOnRead^Read_29/ReadVariableOp^Read_3/DisableCopyOnRead^Read_3/ReadVariableOp^Read_30/DisableCopyOnRead^Read_30/ReadVariableOp^Read_31/DisableCopyOnRead^Read_31/ReadVariableOp^Read_32/DisableCopyOnRead^Read_32/ReadVariableOp^Read_33/DisableCopyOnRead^Read_33/ReadVariableOp^Read_34/DisableCopyOnRead^Read_34/ReadVariableOp^Read_35/DisableCopyOnRead^Read_35/ReadVariableOp^Read_36/DisableCopyOnRead^Read_36/ReadVariableOp^Read_37/DisableCopyOnRead^Read_37/ReadVariableOp^Read_38/DisableCopyOnRead^Read_38/ReadVariableOp^Read_39/DisableCopyOnRead^Read_39/ReadVariableOp^Read_4/DisableCopyOnRead^Read_4/ReadVariableOp^Read_40/DisableCopyOnRead^Read_40/ReadVariableOp^Read_41/DisableCopyOnRead^Read_41/ReadVariableOp^Read_42/DisableCopyOnRead^Read_42/ReadVariableOp^Read_43/DisableCopyOnRead^Read_43/ReadVariableOp^Read_44/DisableCopyOnRead^Read_44/ReadVariableOp^Read_45/DisableCopyOnRead^Read_45/ReadVariableOp^Read_46/DisableCopyOnRead^Read_46/ReadVariableOp^Read_47/DisableCopyOnRead^Read_47/ReadVariableOp^Read_48/DisableCopyOnRead^Read_48/ReadVariableOp^Read_49/DisableCopyOnRead^Read_49/ReadVariableOp^Read_5/DisableCopyOnRead^Read_5/ReadVariableOp^Read_50/DisableCopyOnRead^Read_50/ReadVariableOp^Read_51/DisableCopyOnRead^Read_51/ReadVariableOp^Read_52/DisableCopyOnRead^Read_52/ReadVariableOp^Read_53/DisableCopyOnRead^Read_53/ReadVariableOp^Read_6/DisableCopyOnRead^Read_6/ReadVariableOp^Read_7/DisableCopyOnRead^Read_7/ReadVariableOp^Read_8/DisableCopyOnRead^Read_8/ReadVariableOp^Read_9/DisableCopyOnRead^Read_9/ReadVariableOp*
_output_shapes
 "%
identity_109Identity_109:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesr
p: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2(
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
Read_4/ReadVariableOpRead_4/ReadVariableOp26
Read_40/DisableCopyOnReadRead_40/DisableCopyOnRead20
Read_40/ReadVariableOpRead_40/ReadVariableOp26
Read_41/DisableCopyOnReadRead_41/DisableCopyOnRead20
Read_41/ReadVariableOpRead_41/ReadVariableOp26
Read_42/DisableCopyOnReadRead_42/DisableCopyOnRead20
Read_42/ReadVariableOpRead_42/ReadVariableOp26
Read_43/DisableCopyOnReadRead_43/DisableCopyOnRead20
Read_43/ReadVariableOpRead_43/ReadVariableOp26
Read_44/DisableCopyOnReadRead_44/DisableCopyOnRead20
Read_44/ReadVariableOpRead_44/ReadVariableOp26
Read_45/DisableCopyOnReadRead_45/DisableCopyOnRead20
Read_45/ReadVariableOpRead_45/ReadVariableOp26
Read_46/DisableCopyOnReadRead_46/DisableCopyOnRead20
Read_46/ReadVariableOpRead_46/ReadVariableOp26
Read_47/DisableCopyOnReadRead_47/DisableCopyOnRead20
Read_47/ReadVariableOpRead_47/ReadVariableOp26
Read_48/DisableCopyOnReadRead_48/DisableCopyOnRead20
Read_48/ReadVariableOpRead_48/ReadVariableOp26
Read_49/DisableCopyOnReadRead_49/DisableCopyOnRead20
Read_49/ReadVariableOpRead_49/ReadVariableOp24
Read_5/DisableCopyOnReadRead_5/DisableCopyOnRead2.
Read_5/ReadVariableOpRead_5/ReadVariableOp26
Read_50/DisableCopyOnReadRead_50/DisableCopyOnRead20
Read_50/ReadVariableOpRead_50/ReadVariableOp26
Read_51/DisableCopyOnReadRead_51/DisableCopyOnRead20
Read_51/ReadVariableOpRead_51/ReadVariableOp26
Read_52/DisableCopyOnReadRead_52/DisableCopyOnRead20
Read_52/ReadVariableOpRead_52/ReadVariableOp26
Read_53/DisableCopyOnReadRead_53/DisableCopyOnRead20
Read_53/ReadVariableOpRead_53/ReadVariableOp24
Read_6/DisableCopyOnReadRead_6/DisableCopyOnRead2.
Read_6/ReadVariableOpRead_6/ReadVariableOp24
Read_7/DisableCopyOnReadRead_7/DisableCopyOnRead2.
Read_7/ReadVariableOpRead_7/ReadVariableOp24
Read_8/DisableCopyOnReadRead_8/DisableCopyOnRead2.
Read_8/ReadVariableOpRead_8/ReadVariableOp24
Read_9/DisableCopyOnReadRead_9/DisableCopyOnRead2.
Read_9/ReadVariableOpRead_9/ReadVariableOp:=79

_output_shapes
: 

_user_specified_nameConst:%6!

_user_specified_namecount:%5!

_user_specified_nametotal:'4#
!
_user_specified_name	count_1:'3#
!
_user_specified_name	total_1:D2@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_359/bias:D1@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_359/bias:F0B
@
_user_specified_name(&Adam/v/sparse_model_8/dense_359/kernel:F/B
@
_user_specified_name(&Adam/m/sparse_model_8/dense_359/kernel:R.N
L
_user_specified_name42Adam/v/sparse_model_8/batch_normalization_242/beta:R-N
L
_user_specified_name42Adam/m/sparse_model_8/batch_normalization_242/beta:S,O
M
_user_specified_name53Adam/v/sparse_model_8/batch_normalization_242/gamma:S+O
M
_user_specified_name53Adam/m/sparse_model_8/batch_normalization_242/gamma:D*@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_358/bias:D)@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_358/bias:F(B
@
_user_specified_name(&Adam/v/sparse_model_8/dense_358/kernel:F'B
@
_user_specified_name(&Adam/m/sparse_model_8/dense_358/kernel:R&N
L
_user_specified_name42Adam/v/sparse_model_8/batch_normalization_241/beta:R%N
L
_user_specified_name42Adam/m/sparse_model_8/batch_normalization_241/beta:S$O
M
_user_specified_name53Adam/v/sparse_model_8/batch_normalization_241/gamma:S#O
M
_user_specified_name53Adam/m/sparse_model_8/batch_normalization_241/gamma:D"@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_357/bias:D!@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_357/bias:F B
@
_user_specified_name(&Adam/v/sparse_model_8/dense_357/kernel:FB
@
_user_specified_name(&Adam/m/sparse_model_8/dense_357/kernel:RN
L
_user_specified_name42Adam/v/sparse_model_8/batch_normalization_240/beta:RN
L
_user_specified_name42Adam/m/sparse_model_8/batch_normalization_240/beta:SO
M
_user_specified_name53Adam/v/sparse_model_8/batch_normalization_240/gamma:SO
M
_user_specified_name53Adam/m/sparse_model_8/batch_normalization_240/gamma:D@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_356/bias:D@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_356/bias:FB
@
_user_specified_name(&Adam/v/sparse_model_8/dense_356/kernel:FB
@
_user_specified_name(&Adam/m/sparse_model_8/dense_356/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:=9
7
_user_specified_namesparse_model_8/dense_359/bias:?;
9
_user_specified_name!sparse_model_8/dense_359/kernel:VR
P
_user_specified_name86sparse_model_8/batch_normalization_242/moving_variance:RN
L
_user_specified_name42sparse_model_8/batch_normalization_242/moving_mean:KG
E
_user_specified_name-+sparse_model_8/batch_normalization_242/beta:LH
F
_user_specified_name.,sparse_model_8/batch_normalization_242/gamma:=9
7
_user_specified_namesparse_model_8/dense_358/bias:?;
9
_user_specified_name!sparse_model_8/dense_358/kernel:VR
P
_user_specified_name86sparse_model_8/batch_normalization_241/moving_variance:RN
L
_user_specified_name42sparse_model_8/batch_normalization_241/moving_mean:K
G
E
_user_specified_name-+sparse_model_8/batch_normalization_241/beta:L	H
F
_user_specified_name.,sparse_model_8/batch_normalization_241/gamma:=9
7
_user_specified_namesparse_model_8/dense_357/bias:?;
9
_user_specified_name!sparse_model_8/dense_357/kernel:VR
P
_user_specified_name86sparse_model_8/batch_normalization_240/moving_variance:RN
L
_user_specified_name42sparse_model_8/batch_normalization_240/moving_mean:KG
E
_user_specified_name-+sparse_model_8/batch_normalization_240/beta:LH
F
_user_specified_name.,sparse_model_8/batch_normalization_240/gamma:=9
7
_user_specified_namesparse_model_8/dense_356/bias:?;
9
_user_specified_name!sparse_model_8/dense_356/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
/__inference_sparse_model_8_layer_call_fn_448855
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:

unknown_13:

unknown_14:

unknown_15:

unknown_16:

unknown_17:

unknown_18:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448765o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name448851:&"
 
_user_specified_name448849:&"
 
_user_specified_name448847:&"
 
_user_specified_name448845:&"
 
_user_specified_name448843:&"
 
_user_specified_name448841:&"
 
_user_specified_name448839:&"
 
_user_specified_name448837:&"
 
_user_specified_name448835:&"
 
_user_specified_name448833:&
"
 
_user_specified_name448831:&	"
 
_user_specified_name448829:&"
 
_user_specified_name448827:&"
 
_user_specified_name448825:&"
 
_user_specified_name448823:&"
 
_user_specified_name448821:&"
 
_user_specified_name448819:&"
 
_user_specified_name448817:&"
 
_user_specified_name448815:&"
 
_user_specified_name448813:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�	
�
8__inference_batch_normalization_240_layer_call_fn_448971

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_448416o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name448967:&"
 
_user_specified_name448965:&"
 
_user_specified_name448963:&"
 
_user_specified_name448961:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_448436

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
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
:���������
 
_user_specified_nameinputs
�%
�
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_448416

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
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
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
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
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
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
:���������
 
_user_specified_nameinputs
�
�
/__inference_sparse_model_8_layer_call_fn_448810
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:

unknown_13:

unknown_14:

unknown_15:

unknown_16:

unknown_17:

unknown_18:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448714o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name448806:&"
 
_user_specified_name448804:&"
 
_user_specified_name448802:&"
 
_user_specified_name448800:&"
 
_user_specified_name448798:&"
 
_user_specified_name448796:&"
 
_user_specified_name448794:&"
 
_user_specified_name448792:&"
 
_user_specified_name448790:&"
 
_user_specified_name448788:&
"
 
_user_specified_name448786:&	"
 
_user_specified_name448784:&"
 
_user_specified_name448782:&"
 
_user_specified_name448780:&"
 
_user_specified_name448778:&"
 
_user_specified_name448776:&"
 
_user_specified_name448774:&"
 
_user_specified_name448772:&"
 
_user_specified_name448770:&"
 
_user_specified_name448768:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�%
�
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_448496

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
E__inference_dense_359_layer_call_and_return_conditional_losses_449256

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
�%
�
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449216

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
8__inference_batch_normalization_242_layer_call_fn_449169

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
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_448576o
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
 
_user_specified_name449165:&"
 
_user_specified_name449163:&"
 
_user_specified_name449161:&"
 
_user_specified_name449159:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
*__inference_dense_359_layer_call_fn_449245

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
E__inference_dense_359_layer_call_and_return_conditional_losses_448707o
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
 
_user_specified_name449241:&"
 
_user_specified_name449239:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449137

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
8__inference_batch_normalization_240_layer_call_fn_448984

inputs
unknown:
	unknown_0:
	unknown_1:
	unknown_2:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_448436o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name448980:&"
 
_user_specified_name448978:&"
 
_user_specified_name448976:&"
 
_user_specified_name448974:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�3
�	
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448765
input_1#
dense_356_448717:	�
dense_356_448719:,
batch_normalization_240_448722:,
batch_normalization_240_448724:,
batch_normalization_240_448726:,
batch_normalization_240_448728:"
dense_357_448731:
dense_357_448733:,
batch_normalization_241_448736:,
batch_normalization_241_448738:,
batch_normalization_241_448740:,
batch_normalization_241_448742:"
dense_358_448745:
dense_358_448747:,
batch_normalization_242_448750:,
batch_normalization_242_448752:,
batch_normalization_242_448754:,
batch_normalization_242_448756:"
dense_359_448759:
dense_359_448761:
identity��/batch_normalization_240/StatefulPartitionedCall�/batch_normalization_241/StatefulPartitionedCall�/batch_normalization_242/StatefulPartitionedCall�!dense_356/StatefulPartitionedCall�!dense_357/StatefulPartitionedCall�!dense_358/StatefulPartitionedCall�!dense_359/StatefulPartitionedCall�
!dense_356/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_356_448717dense_356_448719*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_356_layer_call_and_return_conditional_losses_448634�
/batch_normalization_240/StatefulPartitionedCallStatefulPartitionedCall*dense_356/StatefulPartitionedCall:output:0batch_normalization_240_448722batch_normalization_240_448724batch_normalization_240_448726batch_normalization_240_448728*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*&
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_448436�
!dense_357/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_240/StatefulPartitionedCall:output:0dense_357_448731dense_357_448733*
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
E__inference_dense_357_layer_call_and_return_conditional_losses_448658�
/batch_normalization_241/StatefulPartitionedCallStatefulPartitionedCall*dense_357/StatefulPartitionedCall:output:0batch_normalization_241_448736batch_normalization_241_448738batch_normalization_241_448740batch_normalization_241_448742*
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
GPU 2J 8� *\
fWRU
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_448516�
!dense_358/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_241/StatefulPartitionedCall:output:0dense_358_448745dense_358_448747*
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
E__inference_dense_358_layer_call_and_return_conditional_losses_448682�
/batch_normalization_242/StatefulPartitionedCallStatefulPartitionedCall*dense_358/StatefulPartitionedCall:output:0batch_normalization_242_448750batch_normalization_242_448752batch_normalization_242_448754batch_normalization_242_448756*
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
GPU 2J 8� *\
fWRU
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_448596�
!dense_359/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_242/StatefulPartitionedCall:output:0dense_359_448759dense_359_448761*
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
E__inference_dense_359_layer_call_and_return_conditional_losses_448707y
IdentityIdentity*dense_359/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp0^batch_normalization_240/StatefulPartitionedCall0^batch_normalization_241/StatefulPartitionedCall0^batch_normalization_242/StatefulPartitionedCall"^dense_356/StatefulPartitionedCall"^dense_357/StatefulPartitionedCall"^dense_358/StatefulPartitionedCall"^dense_359/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 2b
/batch_normalization_240/StatefulPartitionedCall/batch_normalization_240/StatefulPartitionedCall2b
/batch_normalization_241/StatefulPartitionedCall/batch_normalization_241/StatefulPartitionedCall2b
/batch_normalization_242/StatefulPartitionedCall/batch_normalization_242/StatefulPartitionedCall2F
!dense_356/StatefulPartitionedCall!dense_356/StatefulPartitionedCall2F
!dense_357/StatefulPartitionedCall!dense_357/StatefulPartitionedCall2F
!dense_358/StatefulPartitionedCall!dense_358/StatefulPartitionedCall2F
!dense_359/StatefulPartitionedCall!dense_359/StatefulPartitionedCall:&"
 
_user_specified_name448761:&"
 
_user_specified_name448759:&"
 
_user_specified_name448756:&"
 
_user_specified_name448754:&"
 
_user_specified_name448752:&"
 
_user_specified_name448750:&"
 
_user_specified_name448747:&"
 
_user_specified_name448745:&"
 
_user_specified_name448742:&"
 
_user_specified_name448740:&
"
 
_user_specified_name448738:&	"
 
_user_specified_name448736:&"
 
_user_specified_name448733:&"
 
_user_specified_name448731:&"
 
_user_specified_name448728:&"
 
_user_specified_name448726:&"
 
_user_specified_name448724:&"
 
_user_specified_name448722:&"
 
_user_specified_name448719:&"
 
_user_specified_name448717:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
��
�
!__inference__wrapped_model_448382
input_1J
7sparse_model_8_dense_356_matmul_readvariableop_resource:	�F
8sparse_model_8_dense_356_biasadd_readvariableop_resource:Q
Csparse_model_8_batch_normalization_240_cast_readvariableop_resource:S
Esparse_model_8_batch_normalization_240_cast_1_readvariableop_resource:S
Esparse_model_8_batch_normalization_240_cast_2_readvariableop_resource:S
Esparse_model_8_batch_normalization_240_cast_3_readvariableop_resource:I
7sparse_model_8_dense_357_matmul_readvariableop_resource:F
8sparse_model_8_dense_357_biasadd_readvariableop_resource:Q
Csparse_model_8_batch_normalization_241_cast_readvariableop_resource:S
Esparse_model_8_batch_normalization_241_cast_1_readvariableop_resource:S
Esparse_model_8_batch_normalization_241_cast_2_readvariableop_resource:S
Esparse_model_8_batch_normalization_241_cast_3_readvariableop_resource:I
7sparse_model_8_dense_358_matmul_readvariableop_resource:F
8sparse_model_8_dense_358_biasadd_readvariableop_resource:Q
Csparse_model_8_batch_normalization_242_cast_readvariableop_resource:S
Esparse_model_8_batch_normalization_242_cast_1_readvariableop_resource:S
Esparse_model_8_batch_normalization_242_cast_2_readvariableop_resource:S
Esparse_model_8_batch_normalization_242_cast_3_readvariableop_resource:I
7sparse_model_8_dense_359_matmul_readvariableop_resource:F
8sparse_model_8_dense_359_biasadd_readvariableop_resource:
identity��:sparse_model_8/batch_normalization_240/Cast/ReadVariableOp�<sparse_model_8/batch_normalization_240/Cast_1/ReadVariableOp�<sparse_model_8/batch_normalization_240/Cast_2/ReadVariableOp�<sparse_model_8/batch_normalization_240/Cast_3/ReadVariableOp�:sparse_model_8/batch_normalization_241/Cast/ReadVariableOp�<sparse_model_8/batch_normalization_241/Cast_1/ReadVariableOp�<sparse_model_8/batch_normalization_241/Cast_2/ReadVariableOp�<sparse_model_8/batch_normalization_241/Cast_3/ReadVariableOp�:sparse_model_8/batch_normalization_242/Cast/ReadVariableOp�<sparse_model_8/batch_normalization_242/Cast_1/ReadVariableOp�<sparse_model_8/batch_normalization_242/Cast_2/ReadVariableOp�<sparse_model_8/batch_normalization_242/Cast_3/ReadVariableOp�/sparse_model_8/dense_356/BiasAdd/ReadVariableOp�.sparse_model_8/dense_356/MatMul/ReadVariableOp�/sparse_model_8/dense_357/BiasAdd/ReadVariableOp�.sparse_model_8/dense_357/MatMul/ReadVariableOp�/sparse_model_8/dense_358/BiasAdd/ReadVariableOp�.sparse_model_8/dense_358/MatMul/ReadVariableOp�/sparse_model_8/dense_359/BiasAdd/ReadVariableOp�.sparse_model_8/dense_359/MatMul/ReadVariableOp�
.sparse_model_8/dense_356/MatMul/ReadVariableOpReadVariableOp7sparse_model_8_dense_356_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sparse_model_8/dense_356/MatMulMatMulinput_16sparse_model_8/dense_356/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
/sparse_model_8/dense_356/BiasAdd/ReadVariableOpReadVariableOp8sparse_model_8_dense_356_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 sparse_model_8/dense_356/BiasAddBiasAdd)sparse_model_8/dense_356/MatMul:product:07sparse_model_8/dense_356/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
:sparse_model_8/batch_normalization_240/Cast/ReadVariableOpReadVariableOpCsparse_model_8_batch_normalization_240_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_240/Cast_1/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_240_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_240/Cast_2/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_240_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_240/Cast_3/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_240_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0{
6sparse_model_8/batch_normalization_240/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
4sparse_model_8/batch_normalization_240/batchnorm/addAddV2Dsparse_model_8/batch_normalization_240/Cast_1/ReadVariableOp:value:0?sparse_model_8/batch_normalization_240/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_240/batchnorm/RsqrtRsqrt8sparse_model_8/batch_normalization_240/batchnorm/add:z:0*
T0*
_output_shapes
:�
4sparse_model_8/batch_normalization_240/batchnorm/mulMul:sparse_model_8/batch_normalization_240/batchnorm/Rsqrt:y:0Dsparse_model_8/batch_normalization_240/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_240/batchnorm/mul_1Mul)sparse_model_8/dense_356/BiasAdd:output:08sparse_model_8/batch_normalization_240/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
6sparse_model_8/batch_normalization_240/batchnorm/mul_2MulBsparse_model_8/batch_normalization_240/Cast/ReadVariableOp:value:08sparse_model_8/batch_normalization_240/batchnorm/mul:z:0*
T0*
_output_shapes
:�
4sparse_model_8/batch_normalization_240/batchnorm/subSubDsparse_model_8/batch_normalization_240/Cast_2/ReadVariableOp:value:0:sparse_model_8/batch_normalization_240/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_240/batchnorm/add_1AddV2:sparse_model_8/batch_normalization_240/batchnorm/mul_1:z:08sparse_model_8/batch_normalization_240/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
.sparse_model_8/dense_357/MatMul/ReadVariableOpReadVariableOp7sparse_model_8_dense_357_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
sparse_model_8/dense_357/MatMulMatMul:sparse_model_8/batch_normalization_240/batchnorm/add_1:z:06sparse_model_8/dense_357/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
/sparse_model_8/dense_357/BiasAdd/ReadVariableOpReadVariableOp8sparse_model_8_dense_357_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 sparse_model_8/dense_357/BiasAddBiasAdd)sparse_model_8/dense_357/MatMul:product:07sparse_model_8/dense_357/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
:sparse_model_8/batch_normalization_241/Cast/ReadVariableOpReadVariableOpCsparse_model_8_batch_normalization_241_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_241/Cast_1/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_241_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_241/Cast_2/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_241_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_241/Cast_3/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_241_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0{
6sparse_model_8/batch_normalization_241/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
4sparse_model_8/batch_normalization_241/batchnorm/addAddV2Dsparse_model_8/batch_normalization_241/Cast_1/ReadVariableOp:value:0?sparse_model_8/batch_normalization_241/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_241/batchnorm/RsqrtRsqrt8sparse_model_8/batch_normalization_241/batchnorm/add:z:0*
T0*
_output_shapes
:�
4sparse_model_8/batch_normalization_241/batchnorm/mulMul:sparse_model_8/batch_normalization_241/batchnorm/Rsqrt:y:0Dsparse_model_8/batch_normalization_241/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_241/batchnorm/mul_1Mul)sparse_model_8/dense_357/BiasAdd:output:08sparse_model_8/batch_normalization_241/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
6sparse_model_8/batch_normalization_241/batchnorm/mul_2MulBsparse_model_8/batch_normalization_241/Cast/ReadVariableOp:value:08sparse_model_8/batch_normalization_241/batchnorm/mul:z:0*
T0*
_output_shapes
:�
4sparse_model_8/batch_normalization_241/batchnorm/subSubDsparse_model_8/batch_normalization_241/Cast_2/ReadVariableOp:value:0:sparse_model_8/batch_normalization_241/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_241/batchnorm/add_1AddV2:sparse_model_8/batch_normalization_241/batchnorm/mul_1:z:08sparse_model_8/batch_normalization_241/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
.sparse_model_8/dense_358/MatMul/ReadVariableOpReadVariableOp7sparse_model_8_dense_358_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
sparse_model_8/dense_358/MatMulMatMul:sparse_model_8/batch_normalization_241/batchnorm/add_1:z:06sparse_model_8/dense_358/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
/sparse_model_8/dense_358/BiasAdd/ReadVariableOpReadVariableOp8sparse_model_8_dense_358_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 sparse_model_8/dense_358/BiasAddBiasAdd)sparse_model_8/dense_358/MatMul:product:07sparse_model_8/dense_358/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
:sparse_model_8/batch_normalization_242/Cast/ReadVariableOpReadVariableOpCsparse_model_8_batch_normalization_242_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_242/Cast_1/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_242_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_242/Cast_2/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_242_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
<sparse_model_8/batch_normalization_242/Cast_3/ReadVariableOpReadVariableOpEsparse_model_8_batch_normalization_242_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0{
6sparse_model_8/batch_normalization_242/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
4sparse_model_8/batch_normalization_242/batchnorm/addAddV2Dsparse_model_8/batch_normalization_242/Cast_1/ReadVariableOp:value:0?sparse_model_8/batch_normalization_242/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_242/batchnorm/RsqrtRsqrt8sparse_model_8/batch_normalization_242/batchnorm/add:z:0*
T0*
_output_shapes
:�
4sparse_model_8/batch_normalization_242/batchnorm/mulMul:sparse_model_8/batch_normalization_242/batchnorm/Rsqrt:y:0Dsparse_model_8/batch_normalization_242/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_242/batchnorm/mul_1Mul)sparse_model_8/dense_358/BiasAdd:output:08sparse_model_8/batch_normalization_242/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
6sparse_model_8/batch_normalization_242/batchnorm/mul_2MulBsparse_model_8/batch_normalization_242/Cast/ReadVariableOp:value:08sparse_model_8/batch_normalization_242/batchnorm/mul:z:0*
T0*
_output_shapes
:�
4sparse_model_8/batch_normalization_242/batchnorm/subSubDsparse_model_8/batch_normalization_242/Cast_2/ReadVariableOp:value:0:sparse_model_8/batch_normalization_242/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
6sparse_model_8/batch_normalization_242/batchnorm/add_1AddV2:sparse_model_8/batch_normalization_242/batchnorm/mul_1:z:08sparse_model_8/batch_normalization_242/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
.sparse_model_8/dense_359/MatMul/ReadVariableOpReadVariableOp7sparse_model_8_dense_359_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
sparse_model_8/dense_359/MatMulMatMul:sparse_model_8/batch_normalization_242/batchnorm/add_1:z:06sparse_model_8/dense_359/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
/sparse_model_8/dense_359/BiasAdd/ReadVariableOpReadVariableOp8sparse_model_8_dense_359_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 sparse_model_8/dense_359/BiasAddBiasAdd)sparse_model_8/dense_359/MatMul:product:07sparse_model_8/dense_359/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 sparse_model_8/dense_359/SigmoidSigmoid)sparse_model_8/dense_359/BiasAdd:output:0*
T0*'
_output_shapes
:���������s
IdentityIdentity$sparse_model_8/dense_359/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������	
NoOpNoOp;^sparse_model_8/batch_normalization_240/Cast/ReadVariableOp=^sparse_model_8/batch_normalization_240/Cast_1/ReadVariableOp=^sparse_model_8/batch_normalization_240/Cast_2/ReadVariableOp=^sparse_model_8/batch_normalization_240/Cast_3/ReadVariableOp;^sparse_model_8/batch_normalization_241/Cast/ReadVariableOp=^sparse_model_8/batch_normalization_241/Cast_1/ReadVariableOp=^sparse_model_8/batch_normalization_241/Cast_2/ReadVariableOp=^sparse_model_8/batch_normalization_241/Cast_3/ReadVariableOp;^sparse_model_8/batch_normalization_242/Cast/ReadVariableOp=^sparse_model_8/batch_normalization_242/Cast_1/ReadVariableOp=^sparse_model_8/batch_normalization_242/Cast_2/ReadVariableOp=^sparse_model_8/batch_normalization_242/Cast_3/ReadVariableOp0^sparse_model_8/dense_356/BiasAdd/ReadVariableOp/^sparse_model_8/dense_356/MatMul/ReadVariableOp0^sparse_model_8/dense_357/BiasAdd/ReadVariableOp/^sparse_model_8/dense_357/MatMul/ReadVariableOp0^sparse_model_8/dense_358/BiasAdd/ReadVariableOp/^sparse_model_8/dense_358/MatMul/ReadVariableOp0^sparse_model_8/dense_359/BiasAdd/ReadVariableOp/^sparse_model_8/dense_359/MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 2x
:sparse_model_8/batch_normalization_240/Cast/ReadVariableOp:sparse_model_8/batch_normalization_240/Cast/ReadVariableOp2|
<sparse_model_8/batch_normalization_240/Cast_1/ReadVariableOp<sparse_model_8/batch_normalization_240/Cast_1/ReadVariableOp2|
<sparse_model_8/batch_normalization_240/Cast_2/ReadVariableOp<sparse_model_8/batch_normalization_240/Cast_2/ReadVariableOp2|
<sparse_model_8/batch_normalization_240/Cast_3/ReadVariableOp<sparse_model_8/batch_normalization_240/Cast_3/ReadVariableOp2x
:sparse_model_8/batch_normalization_241/Cast/ReadVariableOp:sparse_model_8/batch_normalization_241/Cast/ReadVariableOp2|
<sparse_model_8/batch_normalization_241/Cast_1/ReadVariableOp<sparse_model_8/batch_normalization_241/Cast_1/ReadVariableOp2|
<sparse_model_8/batch_normalization_241/Cast_2/ReadVariableOp<sparse_model_8/batch_normalization_241/Cast_2/ReadVariableOp2|
<sparse_model_8/batch_normalization_241/Cast_3/ReadVariableOp<sparse_model_8/batch_normalization_241/Cast_3/ReadVariableOp2x
:sparse_model_8/batch_normalization_242/Cast/ReadVariableOp:sparse_model_8/batch_normalization_242/Cast/ReadVariableOp2|
<sparse_model_8/batch_normalization_242/Cast_1/ReadVariableOp<sparse_model_8/batch_normalization_242/Cast_1/ReadVariableOp2|
<sparse_model_8/batch_normalization_242/Cast_2/ReadVariableOp<sparse_model_8/batch_normalization_242/Cast_2/ReadVariableOp2|
<sparse_model_8/batch_normalization_242/Cast_3/ReadVariableOp<sparse_model_8/batch_normalization_242/Cast_3/ReadVariableOp2b
/sparse_model_8/dense_356/BiasAdd/ReadVariableOp/sparse_model_8/dense_356/BiasAdd/ReadVariableOp2`
.sparse_model_8/dense_356/MatMul/ReadVariableOp.sparse_model_8/dense_356/MatMul/ReadVariableOp2b
/sparse_model_8/dense_357/BiasAdd/ReadVariableOp/sparse_model_8/dense_357/BiasAdd/ReadVariableOp2`
.sparse_model_8/dense_357/MatMul/ReadVariableOp.sparse_model_8/dense_357/MatMul/ReadVariableOp2b
/sparse_model_8/dense_358/BiasAdd/ReadVariableOp/sparse_model_8/dense_358/BiasAdd/ReadVariableOp2`
.sparse_model_8/dense_358/MatMul/ReadVariableOp.sparse_model_8/dense_358/MatMul/ReadVariableOp2b
/sparse_model_8/dense_359/BiasAdd/ReadVariableOp/sparse_model_8/dense_359/BiasAdd/ReadVariableOp2`
.sparse_model_8/dense_359/MatMul/ReadVariableOp.sparse_model_8/dense_359/MatMul/ReadVariableOp:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
"
_user_specified_name
resource:($
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
E__inference_dense_357_layer_call_and_return_conditional_losses_448658

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
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
:���������: : 20
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
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449236

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
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_448576

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
E__inference_dense_358_layer_call_and_return_conditional_losses_449156

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
E__inference_dense_356_layer_call_and_return_conditional_losses_448958

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������S
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
�
�
*__inference_dense_358_layer_call_fn_449146

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
E__inference_dense_358_layer_call_and_return_conditional_losses_448682o
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
 
_user_specified_name449142:&"
 
_user_specified_name449140:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�3
�	
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448714
input_1#
dense_356_448635:	�
dense_356_448637:,
batch_normalization_240_448640:,
batch_normalization_240_448642:,
batch_normalization_240_448644:,
batch_normalization_240_448646:"
dense_357_448659:
dense_357_448661:,
batch_normalization_241_448664:,
batch_normalization_241_448666:,
batch_normalization_241_448668:,
batch_normalization_241_448670:"
dense_358_448683:
dense_358_448685:,
batch_normalization_242_448688:,
batch_normalization_242_448690:,
batch_normalization_242_448692:,
batch_normalization_242_448694:"
dense_359_448708:
dense_359_448710:
identity��/batch_normalization_240/StatefulPartitionedCall�/batch_normalization_241/StatefulPartitionedCall�/batch_normalization_242/StatefulPartitionedCall�!dense_356/StatefulPartitionedCall�!dense_357/StatefulPartitionedCall�!dense_358/StatefulPartitionedCall�!dense_359/StatefulPartitionedCall�
!dense_356/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_356_448635dense_356_448637*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_356_layer_call_and_return_conditional_losses_448634�
/batch_normalization_240/StatefulPartitionedCallStatefulPartitionedCall*dense_356/StatefulPartitionedCall:output:0batch_normalization_240_448640batch_normalization_240_448642batch_normalization_240_448644batch_normalization_240_448646*
Tin	
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *\
fWRU
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_448416�
!dense_357/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_240/StatefulPartitionedCall:output:0dense_357_448659dense_357_448661*
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
E__inference_dense_357_layer_call_and_return_conditional_losses_448658�
/batch_normalization_241/StatefulPartitionedCallStatefulPartitionedCall*dense_357/StatefulPartitionedCall:output:0batch_normalization_241_448664batch_normalization_241_448666batch_normalization_241_448668batch_normalization_241_448670*
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
GPU 2J 8� *\
fWRU
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_448496�
!dense_358/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_241/StatefulPartitionedCall:output:0dense_358_448683dense_358_448685*
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
E__inference_dense_358_layer_call_and_return_conditional_losses_448682�
/batch_normalization_242/StatefulPartitionedCallStatefulPartitionedCall*dense_358/StatefulPartitionedCall:output:0batch_normalization_242_448688batch_normalization_242_448690batch_normalization_242_448692batch_normalization_242_448694*
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
GPU 2J 8� *\
fWRU
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_448576�
!dense_359/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_242/StatefulPartitionedCall:output:0dense_359_448708dense_359_448710*
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
E__inference_dense_359_layer_call_and_return_conditional_losses_448707y
IdentityIdentity*dense_359/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp0^batch_normalization_240/StatefulPartitionedCall0^batch_normalization_241/StatefulPartitionedCall0^batch_normalization_242/StatefulPartitionedCall"^dense_356/StatefulPartitionedCall"^dense_357/StatefulPartitionedCall"^dense_358/StatefulPartitionedCall"^dense_359/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 2b
/batch_normalization_240/StatefulPartitionedCall/batch_normalization_240/StatefulPartitionedCall2b
/batch_normalization_241/StatefulPartitionedCall/batch_normalization_241/StatefulPartitionedCall2b
/batch_normalization_242/StatefulPartitionedCall/batch_normalization_242/StatefulPartitionedCall2F
!dense_356/StatefulPartitionedCall!dense_356/StatefulPartitionedCall2F
!dense_357/StatefulPartitionedCall!dense_357/StatefulPartitionedCall2F
!dense_358/StatefulPartitionedCall!dense_358/StatefulPartitionedCall2F
!dense_359/StatefulPartitionedCall!dense_359/StatefulPartitionedCall:&"
 
_user_specified_name448710:&"
 
_user_specified_name448708:&"
 
_user_specified_name448694:&"
 
_user_specified_name448692:&"
 
_user_specified_name448690:&"
 
_user_specified_name448688:&"
 
_user_specified_name448685:&"
 
_user_specified_name448683:&"
 
_user_specified_name448670:&"
 
_user_specified_name448668:&
"
 
_user_specified_name448666:&	"
 
_user_specified_name448664:&"
 
_user_specified_name448661:&"
 
_user_specified_name448659:&"
 
_user_specified_name448646:&"
 
_user_specified_name448644:&"
 
_user_specified_name448642:&"
 
_user_specified_name448640:&"
 
_user_specified_name448637:&"
 
_user_specified_name448635:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
��
�)
"__inference__traced_restore_449773
file_prefixC
0assignvariableop_sparse_model_8_dense_356_kernel:	�>
0assignvariableop_1_sparse_model_8_dense_356_bias:M
?assignvariableop_2_sparse_model_8_batch_normalization_240_gamma:L
>assignvariableop_3_sparse_model_8_batch_normalization_240_beta:S
Eassignvariableop_4_sparse_model_8_batch_normalization_240_moving_mean:W
Iassignvariableop_5_sparse_model_8_batch_normalization_240_moving_variance:D
2assignvariableop_6_sparse_model_8_dense_357_kernel:>
0assignvariableop_7_sparse_model_8_dense_357_bias:M
?assignvariableop_8_sparse_model_8_batch_normalization_241_gamma:L
>assignvariableop_9_sparse_model_8_batch_normalization_241_beta:T
Fassignvariableop_10_sparse_model_8_batch_normalization_241_moving_mean:X
Jassignvariableop_11_sparse_model_8_batch_normalization_241_moving_variance:E
3assignvariableop_12_sparse_model_8_dense_358_kernel:?
1assignvariableop_13_sparse_model_8_dense_358_bias:N
@assignvariableop_14_sparse_model_8_batch_normalization_242_gamma:M
?assignvariableop_15_sparse_model_8_batch_normalization_242_beta:T
Fassignvariableop_16_sparse_model_8_batch_normalization_242_moving_mean:X
Jassignvariableop_17_sparse_model_8_batch_normalization_242_moving_variance:E
3assignvariableop_18_sparse_model_8_dense_359_kernel:?
1assignvariableop_19_sparse_model_8_dense_359_bias:'
assignvariableop_20_iteration:	 +
!assignvariableop_21_learning_rate: M
:assignvariableop_22_adam_m_sparse_model_8_dense_356_kernel:	�M
:assignvariableop_23_adam_v_sparse_model_8_dense_356_kernel:	�F
8assignvariableop_24_adam_m_sparse_model_8_dense_356_bias:F
8assignvariableop_25_adam_v_sparse_model_8_dense_356_bias:U
Gassignvariableop_26_adam_m_sparse_model_8_batch_normalization_240_gamma:U
Gassignvariableop_27_adam_v_sparse_model_8_batch_normalization_240_gamma:T
Fassignvariableop_28_adam_m_sparse_model_8_batch_normalization_240_beta:T
Fassignvariableop_29_adam_v_sparse_model_8_batch_normalization_240_beta:L
:assignvariableop_30_adam_m_sparse_model_8_dense_357_kernel:L
:assignvariableop_31_adam_v_sparse_model_8_dense_357_kernel:F
8assignvariableop_32_adam_m_sparse_model_8_dense_357_bias:F
8assignvariableop_33_adam_v_sparse_model_8_dense_357_bias:U
Gassignvariableop_34_adam_m_sparse_model_8_batch_normalization_241_gamma:U
Gassignvariableop_35_adam_v_sparse_model_8_batch_normalization_241_gamma:T
Fassignvariableop_36_adam_m_sparse_model_8_batch_normalization_241_beta:T
Fassignvariableop_37_adam_v_sparse_model_8_batch_normalization_241_beta:L
:assignvariableop_38_adam_m_sparse_model_8_dense_358_kernel:L
:assignvariableop_39_adam_v_sparse_model_8_dense_358_kernel:F
8assignvariableop_40_adam_m_sparse_model_8_dense_358_bias:F
8assignvariableop_41_adam_v_sparse_model_8_dense_358_bias:U
Gassignvariableop_42_adam_m_sparse_model_8_batch_normalization_242_gamma:U
Gassignvariableop_43_adam_v_sparse_model_8_batch_normalization_242_gamma:T
Fassignvariableop_44_adam_m_sparse_model_8_batch_normalization_242_beta:T
Fassignvariableop_45_adam_v_sparse_model_8_batch_normalization_242_beta:L
:assignvariableop_46_adam_m_sparse_model_8_dense_359_kernel:L
:assignvariableop_47_adam_v_sparse_model_8_dense_359_kernel:F
8assignvariableop_48_adam_m_sparse_model_8_dense_359_bias:F
8assignvariableop_49_adam_v_sparse_model_8_dense_359_bias:%
assignvariableop_50_total_1: %
assignvariableop_51_count_1: #
assignvariableop_52_total: #
assignvariableop_53_count: 
identity_55��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_44�AssignVariableOp_45�AssignVariableOp_46�AssignVariableOp_47�AssignVariableOp_48�AssignVariableOp_49�AssignVariableOp_5�AssignVariableOp_50�AssignVariableOp_51�AssignVariableOp_52�AssignVariableOp_53�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:7*
dtype0*�
value�B�7B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB'variables/18/.ATTRIBUTES/VARIABLE_VALUEB'variables/19/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:7*
dtype0*�
valuexBv7B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�:::::::::::::::::::::::::::::::::::::::::::::::::::::::*E
dtypes;
927	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOp0assignvariableop_sparse_model_8_dense_356_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp0assignvariableop_1_sparse_model_8_dense_356_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOp?assignvariableop_2_sparse_model_8_batch_normalization_240_gammaIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOp>assignvariableop_3_sparse_model_8_batch_normalization_240_betaIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpEassignvariableop_4_sparse_model_8_batch_normalization_240_moving_meanIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpIassignvariableop_5_sparse_model_8_batch_normalization_240_moving_varianceIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp2assignvariableop_6_sparse_model_8_dense_357_kernelIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp0assignvariableop_7_sparse_model_8_dense_357_biasIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOp?assignvariableop_8_sparse_model_8_batch_normalization_241_gammaIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp>assignvariableop_9_sparse_model_8_batch_normalization_241_betaIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOpFassignvariableop_10_sparse_model_8_batch_normalization_241_moving_meanIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpJassignvariableop_11_sparse_model_8_batch_normalization_241_moving_varianceIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp3assignvariableop_12_sparse_model_8_dense_358_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp1assignvariableop_13_sparse_model_8_dense_358_biasIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOp@assignvariableop_14_sparse_model_8_batch_normalization_242_gammaIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOp?assignvariableop_15_sparse_model_8_batch_normalization_242_betaIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOpFassignvariableop_16_sparse_model_8_batch_normalization_242_moving_meanIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOpJassignvariableop_17_sparse_model_8_batch_normalization_242_moving_varianceIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp3assignvariableop_18_sparse_model_8_dense_359_kernelIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp1assignvariableop_19_sparse_model_8_dense_359_biasIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_20AssignVariableOpassignvariableop_20_iterationIdentity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOp!assignvariableop_21_learning_rateIdentity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOp:assignvariableop_22_adam_m_sparse_model_8_dense_356_kernelIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOp:assignvariableop_23_adam_v_sparse_model_8_dense_356_kernelIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOp8assignvariableop_24_adam_m_sparse_model_8_dense_356_biasIdentity_24:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOp8assignvariableop_25_adam_v_sparse_model_8_dense_356_biasIdentity_25:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpGassignvariableop_26_adam_m_sparse_model_8_batch_normalization_240_gammaIdentity_26:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOpGassignvariableop_27_adam_v_sparse_model_8_batch_normalization_240_gammaIdentity_27:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpFassignvariableop_28_adam_m_sparse_model_8_batch_normalization_240_betaIdentity_28:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpFassignvariableop_29_adam_v_sparse_model_8_batch_normalization_240_betaIdentity_29:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOp:assignvariableop_30_adam_m_sparse_model_8_dense_357_kernelIdentity_30:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOp:assignvariableop_31_adam_v_sparse_model_8_dense_357_kernelIdentity_31:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOp8assignvariableop_32_adam_m_sparse_model_8_dense_357_biasIdentity_32:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOp8assignvariableop_33_adam_v_sparse_model_8_dense_357_biasIdentity_33:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOpGassignvariableop_34_adam_m_sparse_model_8_batch_normalization_241_gammaIdentity_34:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOpGassignvariableop_35_adam_v_sparse_model_8_batch_normalization_241_gammaIdentity_35:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_36AssignVariableOpFassignvariableop_36_adam_m_sparse_model_8_batch_normalization_241_betaIdentity_36:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_37AssignVariableOpFassignvariableop_37_adam_v_sparse_model_8_batch_normalization_241_betaIdentity_37:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_38AssignVariableOp:assignvariableop_38_adam_m_sparse_model_8_dense_358_kernelIdentity_38:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_39AssignVariableOp:assignvariableop_39_adam_v_sparse_model_8_dense_358_kernelIdentity_39:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_40AssignVariableOp8assignvariableop_40_adam_m_sparse_model_8_dense_358_biasIdentity_40:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_41AssignVariableOp8assignvariableop_41_adam_v_sparse_model_8_dense_358_biasIdentity_41:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_42AssignVariableOpGassignvariableop_42_adam_m_sparse_model_8_batch_normalization_242_gammaIdentity_42:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_43AssignVariableOpGassignvariableop_43_adam_v_sparse_model_8_batch_normalization_242_gammaIdentity_43:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_44AssignVariableOpFassignvariableop_44_adam_m_sparse_model_8_batch_normalization_242_betaIdentity_44:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_45AssignVariableOpFassignvariableop_45_adam_v_sparse_model_8_batch_normalization_242_betaIdentity_45:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_46AssignVariableOp:assignvariableop_46_adam_m_sparse_model_8_dense_359_kernelIdentity_46:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_47AssignVariableOp:assignvariableop_47_adam_v_sparse_model_8_dense_359_kernelIdentity_47:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_48AssignVariableOp8assignvariableop_48_adam_m_sparse_model_8_dense_359_biasIdentity_48:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_49AssignVariableOp8assignvariableop_49_adam_v_sparse_model_8_dense_359_biasIdentity_49:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_50AssignVariableOpassignvariableop_50_total_1Identity_50:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_51AssignVariableOpassignvariableop_51_count_1Identity_51:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_52AssignVariableOpassignvariableop_52_totalIdentity_52:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_53AssignVariableOpassignvariableop_53_countIdentity_53:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0Y
NoOpNoOp"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 �	
Identity_54Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_55IdentityIdentity_54:output:0^NoOp_1*
T0*
_output_shapes
: �	
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
_output_shapes
 "#
identity_55Identity_55:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesp
n: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2*
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
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92$
AssignVariableOpAssignVariableOp:%6!

_user_specified_namecount:%5!

_user_specified_nametotal:'4#
!
_user_specified_name	count_1:'3#
!
_user_specified_name	total_1:D2@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_359/bias:D1@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_359/bias:F0B
@
_user_specified_name(&Adam/v/sparse_model_8/dense_359/kernel:F/B
@
_user_specified_name(&Adam/m/sparse_model_8/dense_359/kernel:R.N
L
_user_specified_name42Adam/v/sparse_model_8/batch_normalization_242/beta:R-N
L
_user_specified_name42Adam/m/sparse_model_8/batch_normalization_242/beta:S,O
M
_user_specified_name53Adam/v/sparse_model_8/batch_normalization_242/gamma:S+O
M
_user_specified_name53Adam/m/sparse_model_8/batch_normalization_242/gamma:D*@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_358/bias:D)@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_358/bias:F(B
@
_user_specified_name(&Adam/v/sparse_model_8/dense_358/kernel:F'B
@
_user_specified_name(&Adam/m/sparse_model_8/dense_358/kernel:R&N
L
_user_specified_name42Adam/v/sparse_model_8/batch_normalization_241/beta:R%N
L
_user_specified_name42Adam/m/sparse_model_8/batch_normalization_241/beta:S$O
M
_user_specified_name53Adam/v/sparse_model_8/batch_normalization_241/gamma:S#O
M
_user_specified_name53Adam/m/sparse_model_8/batch_normalization_241/gamma:D"@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_357/bias:D!@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_357/bias:F B
@
_user_specified_name(&Adam/v/sparse_model_8/dense_357/kernel:FB
@
_user_specified_name(&Adam/m/sparse_model_8/dense_357/kernel:RN
L
_user_specified_name42Adam/v/sparse_model_8/batch_normalization_240/beta:RN
L
_user_specified_name42Adam/m/sparse_model_8/batch_normalization_240/beta:SO
M
_user_specified_name53Adam/v/sparse_model_8/batch_normalization_240/gamma:SO
M
_user_specified_name53Adam/m/sparse_model_8/batch_normalization_240/gamma:D@
>
_user_specified_name&$Adam/v/sparse_model_8/dense_356/bias:D@
>
_user_specified_name&$Adam/m/sparse_model_8/dense_356/bias:FB
@
_user_specified_name(&Adam/v/sparse_model_8/dense_356/kernel:FB
@
_user_specified_name(&Adam/m/sparse_model_8/dense_356/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:=9
7
_user_specified_namesparse_model_8/dense_359/bias:?;
9
_user_specified_name!sparse_model_8/dense_359/kernel:VR
P
_user_specified_name86sparse_model_8/batch_normalization_242/moving_variance:RN
L
_user_specified_name42sparse_model_8/batch_normalization_242/moving_mean:KG
E
_user_specified_name-+sparse_model_8/batch_normalization_242/beta:LH
F
_user_specified_name.,sparse_model_8/batch_normalization_242/gamma:=9
7
_user_specified_namesparse_model_8/dense_358/bias:?;
9
_user_specified_name!sparse_model_8/dense_358/kernel:VR
P
_user_specified_name86sparse_model_8/batch_normalization_241/moving_variance:RN
L
_user_specified_name42sparse_model_8/batch_normalization_241/moving_mean:K
G
E
_user_specified_name-+sparse_model_8/batch_normalization_241/beta:L	H
F
_user_specified_name.,sparse_model_8/batch_normalization_241/gamma:=9
7
_user_specified_namesparse_model_8/dense_357/bias:?;
9
_user_specified_name!sparse_model_8/dense_357/kernel:VR
P
_user_specified_name86sparse_model_8/batch_normalization_240/moving_variance:RN
L
_user_specified_name42sparse_model_8/batch_normalization_240/moving_mean:KG
E
_user_specified_name-+sparse_model_8/batch_normalization_240/beta:LH
F
_user_specified_name.,sparse_model_8/batch_normalization_240/gamma:=9
7
_user_specified_namesparse_model_8/dense_356/bias:?;
9
_user_specified_name!sparse_model_8/dense_356/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
$__inference_signature_wrapper_448939
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:

unknown_13:

unknown_14:

unknown_15:

unknown_16:

unknown_17:

unknown_18:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference__wrapped_model_448382o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name448935:&"
 
_user_specified_name448933:&"
 
_user_specified_name448931:&"
 
_user_specified_name448929:&"
 
_user_specified_name448927:&"
 
_user_specified_name448925:&"
 
_user_specified_name448923:&"
 
_user_specified_name448921:&"
 
_user_specified_name448919:&"
 
_user_specified_name448917:&
"
 
_user_specified_name448915:&	"
 
_user_specified_name448913:&"
 
_user_specified_name448911:&"
 
_user_specified_name448909:&"
 
_user_specified_name448907:&"
 
_user_specified_name448905:&"
 
_user_specified_name448903:&"
 
_user_specified_name448901:&"
 
_user_specified_name448899:&"
 
_user_specified_name448897:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�%
�
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449018

inputs5
'assignmovingavg_readvariableop_resource:7
)assignmovingavg_1_readvariableop_resource:*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:
identity��AssignMovingAvg�AssignMovingAvg/ReadVariableOp�AssignMovingAvg_1� AssignMovingAvg_1/ReadVariableOp�Cast/ReadVariableOp�Cast_1/ReadVariableOph
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
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
:*
dtype0�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes
:x
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes
:�
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
:*
dtype0�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes
:~
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes
:�
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype0l
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_1/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������h
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:k
batchnorm/subSubCast/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2@
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
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449038

inputs*
cast_readvariableop_resource:,
cast_1_readvariableop_resource:,
cast_2_readvariableop_resource:,
cast_3_readvariableop_resource:
identity��Cast/ReadVariableOp�Cast_1/ReadVariableOp�Cast_2/ReadVariableOp�Cast_3/ReadVariableOpl
Cast/ReadVariableOpReadVariableOpcast_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_1/ReadVariableOpReadVariableOpcast_1_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_2/ReadVariableOpReadVariableOpcast_2_readvariableop_resource*
_output_shapes
:*
dtype0p
Cast_3/ReadVariableOpReadVariableOpcast_3_readvariableop_resource*
_output_shapes
:*
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
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:m
batchnorm/mulMulbatchnorm/Rsqrt:y:0Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������k
batchnorm/mul_2MulCast/ReadVariableOp:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:m
batchnorm/subSubCast_2/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:r
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*'
_output_shapes
:���������b
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^Cast/ReadVariableOp^Cast_1/ReadVariableOp^Cast_2/ReadVariableOp^Cast_3/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������: : : : 2*
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
:���������
 
_user_specified_nameinputs
�	
�
8__inference_batch_normalization_241_layer_call_fn_449083

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
GPU 2J 8� *\
fWRU
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_448516o
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
 
_user_specified_name449079:&"
 
_user_specified_name449077:&"
 
_user_specified_name449075:&"
 
_user_specified_name449073:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
*__inference_dense_357_layer_call_fn_449047

inputs
unknown:
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
E__inference_dense_357_layer_call_and_return_conditional_losses_448658o
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
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name449043:&"
 
_user_specified_name449041:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_448596

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
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449117

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
�
�
*__inference_dense_356_layer_call_fn_448948

inputs
unknown:	�
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *N
fIRG
E__inference_dense_356_layer_call_and_return_conditional_losses_448634o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������<
NoOpNoOp^StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name448944:&"
 
_user_specified_name448942:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�

�
E__inference_dense_359_layer_call_and_return_conditional_losses_448707

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
�	
�
E__inference_dense_357_layer_call_and_return_conditional_losses_449057

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
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
:���������: : 20
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
:���������
 
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


dense1
	norm2

dense2
	norm3
out
	optimizer

signatures"
_tf_keras_model
�
0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
 15
!16
"17
#18
$19"
trackable_list_wrapper
�
0
1
2
3
4
5
6
7
8
9
10
 11
#12
$13"
trackable_list_wrapper
 "
trackable_list_wrapper
�
%non_trainable_variables

&layers
'metrics
(layer_regularization_losses
)layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�
*trace_0
+trace_12�
/__inference_sparse_model_8_layer_call_fn_448810
/__inference_sparse_model_8_layer_call_fn_448855�
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
 z*trace_0z+trace_1
�
,trace_0
-trace_12�
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448714
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448765�
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
 z,trace_0z-trace_1
�B�
!__inference__wrapped_model_448382input_1"�
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
.	variables
/trainable_variables
0regularization_losses
1	keras_api
2__call__
*3&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
4	variables
5trainable_variables
6regularization_losses
7	keras_api
8__call__
*9&call_and_return_all_conditional_losses
:axis
	gamma
beta
moving_mean
moving_variance"
_tf_keras_layer
�
;	variables
<trainable_variables
=regularization_losses
>	keras_api
?__call__
*@&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
A	variables
Btrainable_variables
Cregularization_losses
D	keras_api
E__call__
*F&call_and_return_all_conditional_losses
Gaxis
	gamma
beta
moving_mean
moving_variance"
_tf_keras_layer
�
H	variables
Itrainable_variables
Jregularization_losses
K	keras_api
L__call__
*M&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
�
N	variables
Otrainable_variables
Pregularization_losses
Q	keras_api
R__call__
*S&call_and_return_all_conditional_losses
Taxis
	gamma
 beta
!moving_mean
"moving_variance"
_tf_keras_layer
�
U	variables
Vtrainable_variables
Wregularization_losses
X	keras_api
Y__call__
*Z&call_and_return_all_conditional_losses

#kernel
$bias"
_tf_keras_layer
�
[
_variables
\_iterations
]_learning_rate
^_index_dict
_
_momentums
`_velocities
a_update_step_xla"
experimentalOptimizer
,
bserving_default"
signature_map
2:0	�2sparse_model_8/dense_356/kernel
+:)2sparse_model_8/dense_356/bias
::82,sparse_model_8/batch_normalization_240/gamma
9:72+sparse_model_8/batch_normalization_240/beta
B:@ (22sparse_model_8/batch_normalization_240/moving_mean
F:D (26sparse_model_8/batch_normalization_240/moving_variance
1:/2sparse_model_8/dense_357/kernel
+:)2sparse_model_8/dense_357/bias
::82,sparse_model_8/batch_normalization_241/gamma
9:72+sparse_model_8/batch_normalization_241/beta
B:@ (22sparse_model_8/batch_normalization_241/moving_mean
F:D (26sparse_model_8/batch_normalization_241/moving_variance
1:/2sparse_model_8/dense_358/kernel
+:)2sparse_model_8/dense_358/bias
::82,sparse_model_8/batch_normalization_242/gamma
9:72+sparse_model_8/batch_normalization_242/beta
B:@ (22sparse_model_8/batch_normalization_242/moving_mean
F:D (26sparse_model_8/batch_normalization_242/moving_variance
1:/2sparse_model_8/dense_359/kernel
+:)2sparse_model_8/dense_359/bias
J
0
1
2
3
!4
"5"
trackable_list_wrapper
Q
0
	1

2
3
4
5
6"
trackable_list_wrapper
.
c0
d1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
/__inference_sparse_model_8_layer_call_fn_448810input_1"�
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
/__inference_sparse_model_8_layer_call_fn_448855input_1"�
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
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448714input_1"�
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
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448765input_1"�
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
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
enon_trainable_variables

flayers
gmetrics
hlayer_regularization_losses
ilayer_metrics
.	variables
/trainable_variables
0regularization_losses
2__call__
*3&call_and_return_all_conditional_losses
&3"call_and_return_conditional_losses"
_generic_user_object
�
jtrace_02�
*__inference_dense_356_layer_call_fn_448948�
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
 zjtrace_0
�
ktrace_02�
E__inference_dense_356_layer_call_and_return_conditional_losses_448958�
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
 zktrace_0
<
0
1
2
3"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
lnon_trainable_variables

mlayers
nmetrics
olayer_regularization_losses
player_metrics
4	variables
5trainable_variables
6regularization_losses
8__call__
*9&call_and_return_all_conditional_losses
&9"call_and_return_conditional_losses"
_generic_user_object
�
qtrace_0
rtrace_12�
8__inference_batch_normalization_240_layer_call_fn_448971
8__inference_batch_normalization_240_layer_call_fn_448984�
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
 zqtrace_0zrtrace_1
�
strace_0
ttrace_12�
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449018
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449038�
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
 zstrace_0zttrace_1
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
unon_trainable_variables

vlayers
wmetrics
xlayer_regularization_losses
ylayer_metrics
;	variables
<trainable_variables
=regularization_losses
?__call__
*@&call_and_return_all_conditional_losses
&@"call_and_return_conditional_losses"
_generic_user_object
�
ztrace_02�
*__inference_dense_357_layer_call_fn_449047�
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
 zztrace_0
�
{trace_02�
E__inference_dense_357_layer_call_and_return_conditional_losses_449057�
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
 z{trace_0
<
0
1
2
3"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
|non_trainable_variables

}layers
~metrics
layer_regularization_losses
�layer_metrics
A	variables
Btrainable_variables
Cregularization_losses
E__call__
*F&call_and_return_all_conditional_losses
&F"call_and_return_conditional_losses"
_generic_user_object
�
�trace_0
�trace_12�
8__inference_batch_normalization_241_layer_call_fn_449070
8__inference_batch_normalization_241_layer_call_fn_449083�
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
 z�trace_0z�trace_1
�
�trace_0
�trace_12�
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449117
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449137�
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
 z�trace_0z�trace_1
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
H	variables
Itrainable_variables
Jregularization_losses
L__call__
*M&call_and_return_all_conditional_losses
&M"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
*__inference_dense_358_layer_call_fn_449146�
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
 z�trace_0
�
�trace_02�
E__inference_dense_358_layer_call_and_return_conditional_losses_449156�
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
 z�trace_0
<
0
 1
!2
"3"
trackable_list_wrapper
.
0
 1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
N	variables
Otrainable_variables
Pregularization_losses
R__call__
*S&call_and_return_all_conditional_losses
&S"call_and_return_conditional_losses"
_generic_user_object
�
�trace_0
�trace_12�
8__inference_batch_normalization_242_layer_call_fn_449169
8__inference_batch_normalization_242_layer_call_fn_449182�
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
 z�trace_0z�trace_1
�
�trace_0
�trace_12�
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449216
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449236�
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
 z�trace_0z�trace_1
 "
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
U	variables
Vtrainable_variables
Wregularization_losses
Y__call__
*Z&call_and_return_all_conditional_losses
&Z"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
*__inference_dense_359_layer_call_fn_449245�
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
 z�trace_0
�
�trace_02�
E__inference_dense_359_layer_call_and_return_conditional_losses_449256�
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
 z�trace_0
�
\0
�1
�2
�3
�4
�5
�6
�7
�8
�9
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
�20
�21
�22
�23
�24
�25
�26
�27
�28"
trackable_list_wrapper
:	 2	iteration
: 2learning_rate
 "
trackable_dict_wrapper
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13"
trackable_list_wrapper
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13"
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
$__inference_signature_wrapper_448939input_1"�
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
*__inference_dense_356_layer_call_fn_448948inputs"�
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
E__inference_dense_356_layer_call_and_return_conditional_losses_448958inputs"�
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
0
1"
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
8__inference_batch_normalization_240_layer_call_fn_448971inputs"�
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
8__inference_batch_normalization_240_layer_call_fn_448984inputs"�
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
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449018inputs"�
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
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449038inputs"�
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
*__inference_dense_357_layer_call_fn_449047inputs"�
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
E__inference_dense_357_layer_call_and_return_conditional_losses_449057inputs"�
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
0
1"
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
8__inference_batch_normalization_241_layer_call_fn_449070inputs"�
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
8__inference_batch_normalization_241_layer_call_fn_449083inputs"�
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
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449117inputs"�
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
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449137inputs"�
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
*__inference_dense_358_layer_call_fn_449146inputs"�
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
E__inference_dense_358_layer_call_and_return_conditional_losses_449156inputs"�
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
!0
"1"
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
8__inference_batch_normalization_242_layer_call_fn_449169inputs"�
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
8__inference_batch_normalization_242_layer_call_fn_449182inputs"�
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
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449216inputs"�
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
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449236inputs"�
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
*__inference_dense_359_layer_call_fn_449245inputs"�
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
E__inference_dense_359_layer_call_and_return_conditional_losses_449256inputs"�
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
7:5	�2&Adam/m/sparse_model_8/dense_356/kernel
7:5	�2&Adam/v/sparse_model_8/dense_356/kernel
0:.2$Adam/m/sparse_model_8/dense_356/bias
0:.2$Adam/v/sparse_model_8/dense_356/bias
?:=23Adam/m/sparse_model_8/batch_normalization_240/gamma
?:=23Adam/v/sparse_model_8/batch_normalization_240/gamma
>:<22Adam/m/sparse_model_8/batch_normalization_240/beta
>:<22Adam/v/sparse_model_8/batch_normalization_240/beta
6:42&Adam/m/sparse_model_8/dense_357/kernel
6:42&Adam/v/sparse_model_8/dense_357/kernel
0:.2$Adam/m/sparse_model_8/dense_357/bias
0:.2$Adam/v/sparse_model_8/dense_357/bias
?:=23Adam/m/sparse_model_8/batch_normalization_241/gamma
?:=23Adam/v/sparse_model_8/batch_normalization_241/gamma
>:<22Adam/m/sparse_model_8/batch_normalization_241/beta
>:<22Adam/v/sparse_model_8/batch_normalization_241/beta
6:42&Adam/m/sparse_model_8/dense_358/kernel
6:42&Adam/v/sparse_model_8/dense_358/kernel
0:.2$Adam/m/sparse_model_8/dense_358/bias
0:.2$Adam/v/sparse_model_8/dense_358/bias
?:=23Adam/m/sparse_model_8/batch_normalization_242/gamma
?:=23Adam/v/sparse_model_8/batch_normalization_242/gamma
>:<22Adam/m/sparse_model_8/batch_normalization_242/beta
>:<22Adam/v/sparse_model_8/batch_normalization_242/beta
6:42&Adam/m/sparse_model_8/dense_359/kernel
6:42&Adam/v/sparse_model_8/dense_359/kernel
0:.2$Adam/m/sparse_model_8/dense_359/bias
0:.2$Adam/v/sparse_model_8/dense_359/bias
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
!__inference__wrapped_model_448382~!" #$1�.
'�$
"�
input_1����������
� "3�0
.
output_1"�
output_1����������
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449018m7�4
-�*
 �
inputs���������
p

 
� ",�)
"�
tensor_0���������
� �
S__inference_batch_normalization_240_layer_call_and_return_conditional_losses_449038m7�4
-�*
 �
inputs���������
p 

 
� ",�)
"�
tensor_0���������
� �
8__inference_batch_normalization_240_layer_call_fn_448971b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_240_layer_call_fn_448984b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449117m7�4
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
S__inference_batch_normalization_241_layer_call_and_return_conditional_losses_449137m7�4
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
8__inference_batch_normalization_241_layer_call_fn_449070b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_241_layer_call_fn_449083b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449216m!" 7�4
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
S__inference_batch_normalization_242_layer_call_and_return_conditional_losses_449236m!" 7�4
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
8__inference_batch_normalization_242_layer_call_fn_449169b!" 7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_242_layer_call_fn_449182b!" 7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
E__inference_dense_356_layer_call_and_return_conditional_losses_448958d0�-
&�#
!�
inputs����������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_356_layer_call_fn_448948Y0�-
&�#
!�
inputs����������
� "!�
unknown����������
E__inference_dense_357_layer_call_and_return_conditional_losses_449057c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_357_layer_call_fn_449047X/�,
%�"
 �
inputs���������
� "!�
unknown����������
E__inference_dense_358_layer_call_and_return_conditional_losses_449156c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_358_layer_call_fn_449146X/�,
%�"
 �
inputs���������
� "!�
unknown����������
E__inference_dense_359_layer_call_and_return_conditional_losses_449256c#$/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_359_layer_call_fn_449245X#$/�,
%�"
 �
inputs���������
� "!�
unknown����������
$__inference_signature_wrapper_448939�!" #$<�9
� 
2�/
-
input_1"�
input_1����������"3�0
.
output_1"�
output_1����������
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448714�!" #$A�>
'�$
"�
input_1����������
�

trainingp",�)
"�
tensor_0���������
� �
J__inference_sparse_model_8_layer_call_and_return_conditional_losses_448765�!" #$A�>
'�$
"�
input_1����������
�

trainingp ",�)
"�
tensor_0���������
� �
/__inference_sparse_model_8_layer_call_fn_448810|!" #$A�>
'�$
"�
input_1����������
�

trainingp"!�
unknown����������
/__inference_sparse_model_8_layer_call_fn_448855|!" #$A�>
'�$
"�
input_1����������
�

trainingp "!�
unknown���������