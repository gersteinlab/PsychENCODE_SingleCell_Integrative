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
 �"serve*2.13.02v2.13.0-rc2-7-g1cb1a030a628�
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
%Adam/v/sparse_model_24/dense_423/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/v/sparse_model_24/dense_423/bias
�
9Adam/v/sparse_model_24/dense_423/bias/Read/ReadVariableOpReadVariableOp%Adam/v/sparse_model_24/dense_423/bias*
_output_shapes
:*
dtype0
�
%Adam/m/sparse_model_24/dense_423/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/m/sparse_model_24/dense_423/bias
�
9Adam/m/sparse_model_24/dense_423/bias/Read/ReadVariableOpReadVariableOp%Adam/m/sparse_model_24/dense_423/bias*
_output_shapes
:*
dtype0
�
'Adam/v/sparse_model_24/dense_423/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/v/sparse_model_24/dense_423/kernel
�
;Adam/v/sparse_model_24/dense_423/kernel/Read/ReadVariableOpReadVariableOp'Adam/v/sparse_model_24/dense_423/kernel*
_output_shapes

:*
dtype0
�
'Adam/m/sparse_model_24/dense_423/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/m/sparse_model_24/dense_423/kernel
�
;Adam/m/sparse_model_24/dense_423/kernel/Read/ReadVariableOpReadVariableOp'Adam/m/sparse_model_24/dense_423/kernel*
_output_shapes

:*
dtype0
�
3Adam/v/sparse_model_24/batch_normalization_290/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/sparse_model_24/batch_normalization_290/beta
�
GAdam/v/sparse_model_24/batch_normalization_290/beta/Read/ReadVariableOpReadVariableOp3Adam/v/sparse_model_24/batch_normalization_290/beta*
_output_shapes
:*
dtype0
�
3Adam/m/sparse_model_24/batch_normalization_290/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/sparse_model_24/batch_normalization_290/beta
�
GAdam/m/sparse_model_24/batch_normalization_290/beta/Read/ReadVariableOpReadVariableOp3Adam/m/sparse_model_24/batch_normalization_290/beta*
_output_shapes
:*
dtype0
�
4Adam/v/sparse_model_24/batch_normalization_290/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/v/sparse_model_24/batch_normalization_290/gamma
�
HAdam/v/sparse_model_24/batch_normalization_290/gamma/Read/ReadVariableOpReadVariableOp4Adam/v/sparse_model_24/batch_normalization_290/gamma*
_output_shapes
:*
dtype0
�
4Adam/m/sparse_model_24/batch_normalization_290/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/m/sparse_model_24/batch_normalization_290/gamma
�
HAdam/m/sparse_model_24/batch_normalization_290/gamma/Read/ReadVariableOpReadVariableOp4Adam/m/sparse_model_24/batch_normalization_290/gamma*
_output_shapes
:*
dtype0
�
%Adam/v/sparse_model_24/dense_422/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/v/sparse_model_24/dense_422/bias
�
9Adam/v/sparse_model_24/dense_422/bias/Read/ReadVariableOpReadVariableOp%Adam/v/sparse_model_24/dense_422/bias*
_output_shapes
:*
dtype0
�
%Adam/m/sparse_model_24/dense_422/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/m/sparse_model_24/dense_422/bias
�
9Adam/m/sparse_model_24/dense_422/bias/Read/ReadVariableOpReadVariableOp%Adam/m/sparse_model_24/dense_422/bias*
_output_shapes
:*
dtype0
�
'Adam/v/sparse_model_24/dense_422/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/v/sparse_model_24/dense_422/kernel
�
;Adam/v/sparse_model_24/dense_422/kernel/Read/ReadVariableOpReadVariableOp'Adam/v/sparse_model_24/dense_422/kernel*
_output_shapes

:*
dtype0
�
'Adam/m/sparse_model_24/dense_422/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/m/sparse_model_24/dense_422/kernel
�
;Adam/m/sparse_model_24/dense_422/kernel/Read/ReadVariableOpReadVariableOp'Adam/m/sparse_model_24/dense_422/kernel*
_output_shapes

:*
dtype0
�
3Adam/v/sparse_model_24/batch_normalization_289/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/sparse_model_24/batch_normalization_289/beta
�
GAdam/v/sparse_model_24/batch_normalization_289/beta/Read/ReadVariableOpReadVariableOp3Adam/v/sparse_model_24/batch_normalization_289/beta*
_output_shapes
:*
dtype0
�
3Adam/m/sparse_model_24/batch_normalization_289/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/sparse_model_24/batch_normalization_289/beta
�
GAdam/m/sparse_model_24/batch_normalization_289/beta/Read/ReadVariableOpReadVariableOp3Adam/m/sparse_model_24/batch_normalization_289/beta*
_output_shapes
:*
dtype0
�
4Adam/v/sparse_model_24/batch_normalization_289/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/v/sparse_model_24/batch_normalization_289/gamma
�
HAdam/v/sparse_model_24/batch_normalization_289/gamma/Read/ReadVariableOpReadVariableOp4Adam/v/sparse_model_24/batch_normalization_289/gamma*
_output_shapes
:*
dtype0
�
4Adam/m/sparse_model_24/batch_normalization_289/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/m/sparse_model_24/batch_normalization_289/gamma
�
HAdam/m/sparse_model_24/batch_normalization_289/gamma/Read/ReadVariableOpReadVariableOp4Adam/m/sparse_model_24/batch_normalization_289/gamma*
_output_shapes
:*
dtype0
�
%Adam/v/sparse_model_24/dense_421/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/v/sparse_model_24/dense_421/bias
�
9Adam/v/sparse_model_24/dense_421/bias/Read/ReadVariableOpReadVariableOp%Adam/v/sparse_model_24/dense_421/bias*
_output_shapes
:*
dtype0
�
%Adam/m/sparse_model_24/dense_421/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/m/sparse_model_24/dense_421/bias
�
9Adam/m/sparse_model_24/dense_421/bias/Read/ReadVariableOpReadVariableOp%Adam/m/sparse_model_24/dense_421/bias*
_output_shapes
:*
dtype0
�
'Adam/v/sparse_model_24/dense_421/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/v/sparse_model_24/dense_421/kernel
�
;Adam/v/sparse_model_24/dense_421/kernel/Read/ReadVariableOpReadVariableOp'Adam/v/sparse_model_24/dense_421/kernel*
_output_shapes

:*
dtype0
�
'Adam/m/sparse_model_24/dense_421/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*8
shared_name)'Adam/m/sparse_model_24/dense_421/kernel
�
;Adam/m/sparse_model_24/dense_421/kernel/Read/ReadVariableOpReadVariableOp'Adam/m/sparse_model_24/dense_421/kernel*
_output_shapes

:*
dtype0
�
3Adam/v/sparse_model_24/batch_normalization_288/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/v/sparse_model_24/batch_normalization_288/beta
�
GAdam/v/sparse_model_24/batch_normalization_288/beta/Read/ReadVariableOpReadVariableOp3Adam/v/sparse_model_24/batch_normalization_288/beta*
_output_shapes
:*
dtype0
�
3Adam/m/sparse_model_24/batch_normalization_288/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53Adam/m/sparse_model_24/batch_normalization_288/beta
�
GAdam/m/sparse_model_24/batch_normalization_288/beta/Read/ReadVariableOpReadVariableOp3Adam/m/sparse_model_24/batch_normalization_288/beta*
_output_shapes
:*
dtype0
�
4Adam/v/sparse_model_24/batch_normalization_288/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/v/sparse_model_24/batch_normalization_288/gamma
�
HAdam/v/sparse_model_24/batch_normalization_288/gamma/Read/ReadVariableOpReadVariableOp4Adam/v/sparse_model_24/batch_normalization_288/gamma*
_output_shapes
:*
dtype0
�
4Adam/m/sparse_model_24/batch_normalization_288/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*E
shared_name64Adam/m/sparse_model_24/batch_normalization_288/gamma
�
HAdam/m/sparse_model_24/batch_normalization_288/gamma/Read/ReadVariableOpReadVariableOp4Adam/m/sparse_model_24/batch_normalization_288/gamma*
_output_shapes
:*
dtype0
�
%Adam/v/sparse_model_24/dense_420/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/v/sparse_model_24/dense_420/bias
�
9Adam/v/sparse_model_24/dense_420/bias/Read/ReadVariableOpReadVariableOp%Adam/v/sparse_model_24/dense_420/bias*
_output_shapes
:*
dtype0
�
%Adam/m/sparse_model_24/dense_420/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%Adam/m/sparse_model_24/dense_420/bias
�
9Adam/m/sparse_model_24/dense_420/bias/Read/ReadVariableOpReadVariableOp%Adam/m/sparse_model_24/dense_420/bias*
_output_shapes
:*
dtype0
�
'Adam/v/sparse_model_24/dense_420/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*8
shared_name)'Adam/v/sparse_model_24/dense_420/kernel
�
;Adam/v/sparse_model_24/dense_420/kernel/Read/ReadVariableOpReadVariableOp'Adam/v/sparse_model_24/dense_420/kernel*
_output_shapes
:	�*
dtype0
�
'Adam/m/sparse_model_24/dense_420/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*8
shared_name)'Adam/m/sparse_model_24/dense_420/kernel
�
;Adam/m/sparse_model_24/dense_420/kernel/Read/ReadVariableOpReadVariableOp'Adam/m/sparse_model_24/dense_420/kernel*
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
sparse_model_24/dense_423/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name sparse_model_24/dense_423/bias
�
2sparse_model_24/dense_423/bias/Read/ReadVariableOpReadVariableOpsparse_model_24/dense_423/bias*
_output_shapes
:*
dtype0
�
 sparse_model_24/dense_423/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*1
shared_name" sparse_model_24/dense_423/kernel
�
4sparse_model_24/dense_423/kernel/Read/ReadVariableOpReadVariableOp sparse_model_24/dense_423/kernel*
_output_shapes

:*
dtype0
�
7sparse_model_24/batch_normalization_290/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*H
shared_name97sparse_model_24/batch_normalization_290/moving_variance
�
Ksparse_model_24/batch_normalization_290/moving_variance/Read/ReadVariableOpReadVariableOp7sparse_model_24/batch_normalization_290/moving_variance*
_output_shapes
:*
dtype0
�
3sparse_model_24/batch_normalization_290/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53sparse_model_24/batch_normalization_290/moving_mean
�
Gsparse_model_24/batch_normalization_290/moving_mean/Read/ReadVariableOpReadVariableOp3sparse_model_24/batch_normalization_290/moving_mean*
_output_shapes
:*
dtype0
�
,sparse_model_24/batch_normalization_290/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,sparse_model_24/batch_normalization_290/beta
�
@sparse_model_24/batch_normalization_290/beta/Read/ReadVariableOpReadVariableOp,sparse_model_24/batch_normalization_290/beta*
_output_shapes
:*
dtype0
�
-sparse_model_24/batch_normalization_290/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*>
shared_name/-sparse_model_24/batch_normalization_290/gamma
�
Asparse_model_24/batch_normalization_290/gamma/Read/ReadVariableOpReadVariableOp-sparse_model_24/batch_normalization_290/gamma*
_output_shapes
:*
dtype0
�
sparse_model_24/dense_422/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name sparse_model_24/dense_422/bias
�
2sparse_model_24/dense_422/bias/Read/ReadVariableOpReadVariableOpsparse_model_24/dense_422/bias*
_output_shapes
:*
dtype0
�
 sparse_model_24/dense_422/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*1
shared_name" sparse_model_24/dense_422/kernel
�
4sparse_model_24/dense_422/kernel/Read/ReadVariableOpReadVariableOp sparse_model_24/dense_422/kernel*
_output_shapes

:*
dtype0
�
7sparse_model_24/batch_normalization_289/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*H
shared_name97sparse_model_24/batch_normalization_289/moving_variance
�
Ksparse_model_24/batch_normalization_289/moving_variance/Read/ReadVariableOpReadVariableOp7sparse_model_24/batch_normalization_289/moving_variance*
_output_shapes
:*
dtype0
�
3sparse_model_24/batch_normalization_289/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53sparse_model_24/batch_normalization_289/moving_mean
�
Gsparse_model_24/batch_normalization_289/moving_mean/Read/ReadVariableOpReadVariableOp3sparse_model_24/batch_normalization_289/moving_mean*
_output_shapes
:*
dtype0
�
,sparse_model_24/batch_normalization_289/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,sparse_model_24/batch_normalization_289/beta
�
@sparse_model_24/batch_normalization_289/beta/Read/ReadVariableOpReadVariableOp,sparse_model_24/batch_normalization_289/beta*
_output_shapes
:*
dtype0
�
-sparse_model_24/batch_normalization_289/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*>
shared_name/-sparse_model_24/batch_normalization_289/gamma
�
Asparse_model_24/batch_normalization_289/gamma/Read/ReadVariableOpReadVariableOp-sparse_model_24/batch_normalization_289/gamma*
_output_shapes
:*
dtype0
�
sparse_model_24/dense_421/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name sparse_model_24/dense_421/bias
�
2sparse_model_24/dense_421/bias/Read/ReadVariableOpReadVariableOpsparse_model_24/dense_421/bias*
_output_shapes
:*
dtype0
�
 sparse_model_24/dense_421/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*1
shared_name" sparse_model_24/dense_421/kernel
�
4sparse_model_24/dense_421/kernel/Read/ReadVariableOpReadVariableOp sparse_model_24/dense_421/kernel*
_output_shapes

:*
dtype0
�
7sparse_model_24/batch_normalization_288/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*H
shared_name97sparse_model_24/batch_normalization_288/moving_variance
�
Ksparse_model_24/batch_normalization_288/moving_variance/Read/ReadVariableOpReadVariableOp7sparse_model_24/batch_normalization_288/moving_variance*
_output_shapes
:*
dtype0
�
3sparse_model_24/batch_normalization_288/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*D
shared_name53sparse_model_24/batch_normalization_288/moving_mean
�
Gsparse_model_24/batch_normalization_288/moving_mean/Read/ReadVariableOpReadVariableOp3sparse_model_24/batch_normalization_288/moving_mean*
_output_shapes
:*
dtype0
�
,sparse_model_24/batch_normalization_288/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*=
shared_name.,sparse_model_24/batch_normalization_288/beta
�
@sparse_model_24/batch_normalization_288/beta/Read/ReadVariableOpReadVariableOp,sparse_model_24/batch_normalization_288/beta*
_output_shapes
:*
dtype0
�
-sparse_model_24/batch_normalization_288/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*>
shared_name/-sparse_model_24/batch_normalization_288/gamma
�
Asparse_model_24/batch_normalization_288/gamma/Read/ReadVariableOpReadVariableOp-sparse_model_24/batch_normalization_288/gamma*
_output_shapes
:*
dtype0
�
sparse_model_24/dense_420/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*/
shared_name sparse_model_24/dense_420/bias
�
2sparse_model_24/dense_420/bias/Read/ReadVariableOpReadVariableOpsparse_model_24/dense_420/bias*
_output_shapes
:*
dtype0
�
 sparse_model_24/dense_420/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*1
shared_name" sparse_model_24/dense_420/kernel
�
4sparse_model_24/dense_420/kernel/Read/ReadVariableOpReadVariableOp sparse_model_24/dense_420/kernel*
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
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1 sparse_model_24/dense_420/kernelsparse_model_24/dense_420/bias3sparse_model_24/batch_normalization_288/moving_mean7sparse_model_24/batch_normalization_288/moving_variance,sparse_model_24/batch_normalization_288/beta-sparse_model_24/batch_normalization_288/gamma sparse_model_24/dense_421/kernelsparse_model_24/dense_421/bias3sparse_model_24/batch_normalization_289/moving_mean7sparse_model_24/batch_normalization_289/moving_variance,sparse_model_24/batch_normalization_289/beta-sparse_model_24/batch_normalization_289/gamma sparse_model_24/dense_422/kernelsparse_model_24/dense_422/bias3sparse_model_24/batch_normalization_290/moving_mean7sparse_model_24/batch_normalization_290/moving_variance,sparse_model_24/batch_normalization_290/beta-sparse_model_24/batch_normalization_290/gamma sparse_model_24/dense_423/kernelsparse_model_24/dense_423/bias* 
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
$__inference_signature_wrapper_525937

NoOpNoOp
�`
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�`
value�`B�` B�`
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
`Z
VARIABLE_VALUE sparse_model_24/dense_420/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEsparse_model_24/dense_420/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
mg
VARIABLE_VALUE-sparse_model_24/batch_normalization_288/gamma&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE,sparse_model_24/batch_normalization_288/beta&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE3sparse_model_24/batch_normalization_288/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
wq
VARIABLE_VALUE7sparse_model_24/batch_normalization_288/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUE sparse_model_24/dense_421/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEsparse_model_24/dense_421/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE*
mg
VARIABLE_VALUE-sparse_model_24/batch_normalization_289/gamma&variables/8/.ATTRIBUTES/VARIABLE_VALUE*
lf
VARIABLE_VALUE,sparse_model_24/batch_normalization_289/beta&variables/9/.ATTRIBUTES/VARIABLE_VALUE*
tn
VARIABLE_VALUE3sparse_model_24/batch_normalization_289/moving_mean'variables/10/.ATTRIBUTES/VARIABLE_VALUE*
xr
VARIABLE_VALUE7sparse_model_24/batch_normalization_289/moving_variance'variables/11/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUE sparse_model_24/dense_422/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEsparse_model_24/dense_422/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE*
nh
VARIABLE_VALUE-sparse_model_24/batch_normalization_290/gamma'variables/14/.ATTRIBUTES/VARIABLE_VALUE*
mg
VARIABLE_VALUE,sparse_model_24/batch_normalization_290/beta'variables/15/.ATTRIBUTES/VARIABLE_VALUE*
tn
VARIABLE_VALUE3sparse_model_24/batch_normalization_290/moving_mean'variables/16/.ATTRIBUTES/VARIABLE_VALUE*
xr
VARIABLE_VALUE7sparse_model_24/batch_normalization_290/moving_variance'variables/17/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUE sparse_model_24/dense_423/kernel'variables/18/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEsparse_model_24/dense_423/bias'variables/19/.ATTRIBUTES/VARIABLE_VALUE*
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
rl
VARIABLE_VALUE'Adam/m/sparse_model_24/dense_420/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE'Adam/v/sparse_model_24/dense_420/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE%Adam/m/sparse_model_24/dense_420/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
pj
VARIABLE_VALUE%Adam/v/sparse_model_24/dense_420/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE4Adam/m/sparse_model_24/batch_normalization_288/gamma1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE4Adam/v/sparse_model_24/batch_normalization_288/gamma1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE3Adam/m/sparse_model_24/batch_normalization_288/beta1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
~x
VARIABLE_VALUE3Adam/v/sparse_model_24/batch_normalization_288/beta1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
rl
VARIABLE_VALUE'Adam/m/sparse_model_24/dense_421/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/v/sparse_model_24/dense_421/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/m/sparse_model_24/dense_421/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/v/sparse_model_24/dense_421/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUE4Adam/m/sparse_model_24/batch_normalization_289/gamma2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUE4Adam/v/sparse_model_24/batch_normalization_289/gamma2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/m/sparse_model_24/batch_normalization_289/beta2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/v/sparse_model_24/batch_normalization_289/beta2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/m/sparse_model_24/dense_422/kernel2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/v/sparse_model_24/dense_422/kernel2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/m/sparse_model_24/dense_422/bias2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/v/sparse_model_24/dense_422/bias2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUE4Adam/m/sparse_model_24/batch_normalization_290/gamma2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUE4Adam/v/sparse_model_24/batch_normalization_290/gamma2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/m/sparse_model_24/batch_normalization_290/beta2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUE3Adam/v/sparse_model_24/batch_normalization_290/beta2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/m/sparse_model_24/dense_423/kernel2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUE*
sm
VARIABLE_VALUE'Adam/v/sparse_model_24/dense_423/kernel2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/m/sparse_model_24/dense_423/bias2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUE*
qk
VARIABLE_VALUE%Adam/v/sparse_model_24/dense_423/bias2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUE*
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
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename sparse_model_24/dense_420/kernelsparse_model_24/dense_420/bias-sparse_model_24/batch_normalization_288/gamma,sparse_model_24/batch_normalization_288/beta3sparse_model_24/batch_normalization_288/moving_mean7sparse_model_24/batch_normalization_288/moving_variance sparse_model_24/dense_421/kernelsparse_model_24/dense_421/bias-sparse_model_24/batch_normalization_289/gamma,sparse_model_24/batch_normalization_289/beta3sparse_model_24/batch_normalization_289/moving_mean7sparse_model_24/batch_normalization_289/moving_variance sparse_model_24/dense_422/kernelsparse_model_24/dense_422/bias-sparse_model_24/batch_normalization_290/gamma,sparse_model_24/batch_normalization_290/beta3sparse_model_24/batch_normalization_290/moving_mean7sparse_model_24/batch_normalization_290/moving_variance sparse_model_24/dense_423/kernelsparse_model_24/dense_423/bias	iterationlearning_rate'Adam/m/sparse_model_24/dense_420/kernel'Adam/v/sparse_model_24/dense_420/kernel%Adam/m/sparse_model_24/dense_420/bias%Adam/v/sparse_model_24/dense_420/bias4Adam/m/sparse_model_24/batch_normalization_288/gamma4Adam/v/sparse_model_24/batch_normalization_288/gamma3Adam/m/sparse_model_24/batch_normalization_288/beta3Adam/v/sparse_model_24/batch_normalization_288/beta'Adam/m/sparse_model_24/dense_421/kernel'Adam/v/sparse_model_24/dense_421/kernel%Adam/m/sparse_model_24/dense_421/bias%Adam/v/sparse_model_24/dense_421/bias4Adam/m/sparse_model_24/batch_normalization_289/gamma4Adam/v/sparse_model_24/batch_normalization_289/gamma3Adam/m/sparse_model_24/batch_normalization_289/beta3Adam/v/sparse_model_24/batch_normalization_289/beta'Adam/m/sparse_model_24/dense_422/kernel'Adam/v/sparse_model_24/dense_422/kernel%Adam/m/sparse_model_24/dense_422/bias%Adam/v/sparse_model_24/dense_422/bias4Adam/m/sparse_model_24/batch_normalization_290/gamma4Adam/v/sparse_model_24/batch_normalization_290/gamma3Adam/m/sparse_model_24/batch_normalization_290/beta3Adam/v/sparse_model_24/batch_normalization_290/beta'Adam/m/sparse_model_24/dense_423/kernel'Adam/v/sparse_model_24/dense_423/kernel%Adam/m/sparse_model_24/dense_423/bias%Adam/v/sparse_model_24/dense_423/biastotal_1count_1totalcountConst*C
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
__inference__traced_save_526600
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename sparse_model_24/dense_420/kernelsparse_model_24/dense_420/bias-sparse_model_24/batch_normalization_288/gamma,sparse_model_24/batch_normalization_288/beta3sparse_model_24/batch_normalization_288/moving_mean7sparse_model_24/batch_normalization_288/moving_variance sparse_model_24/dense_421/kernelsparse_model_24/dense_421/bias-sparse_model_24/batch_normalization_289/gamma,sparse_model_24/batch_normalization_289/beta3sparse_model_24/batch_normalization_289/moving_mean7sparse_model_24/batch_normalization_289/moving_variance sparse_model_24/dense_422/kernelsparse_model_24/dense_422/bias-sparse_model_24/batch_normalization_290/gamma,sparse_model_24/batch_normalization_290/beta3sparse_model_24/batch_normalization_290/moving_mean7sparse_model_24/batch_normalization_290/moving_variance sparse_model_24/dense_423/kernelsparse_model_24/dense_423/bias	iterationlearning_rate'Adam/m/sparse_model_24/dense_420/kernel'Adam/v/sparse_model_24/dense_420/kernel%Adam/m/sparse_model_24/dense_420/bias%Adam/v/sparse_model_24/dense_420/bias4Adam/m/sparse_model_24/batch_normalization_288/gamma4Adam/v/sparse_model_24/batch_normalization_288/gamma3Adam/m/sparse_model_24/batch_normalization_288/beta3Adam/v/sparse_model_24/batch_normalization_288/beta'Adam/m/sparse_model_24/dense_421/kernel'Adam/v/sparse_model_24/dense_421/kernel%Adam/m/sparse_model_24/dense_421/bias%Adam/v/sparse_model_24/dense_421/bias4Adam/m/sparse_model_24/batch_normalization_289/gamma4Adam/v/sparse_model_24/batch_normalization_289/gamma3Adam/m/sparse_model_24/batch_normalization_289/beta3Adam/v/sparse_model_24/batch_normalization_289/beta'Adam/m/sparse_model_24/dense_422/kernel'Adam/v/sparse_model_24/dense_422/kernel%Adam/m/sparse_model_24/dense_422/bias%Adam/v/sparse_model_24/dense_422/bias4Adam/m/sparse_model_24/batch_normalization_290/gamma4Adam/v/sparse_model_24/batch_normalization_290/gamma3Adam/m/sparse_model_24/batch_normalization_290/beta3Adam/v/sparse_model_24/batch_normalization_290/beta'Adam/m/sparse_model_24/dense_423/kernel'Adam/v/sparse_model_24/dense_423/kernel%Adam/m/sparse_model_24/dense_423/bias%Adam/v/sparse_model_24/dense_423/biastotal_1count_1totalcount*B
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
"__inference__traced_restore_526771��
�%
�
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526115

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
�3
�	
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525712
input_1#
dense_420_525633:	�
dense_420_525635:,
batch_normalization_288_525638:,
batch_normalization_288_525640:,
batch_normalization_288_525642:,
batch_normalization_288_525644:"
dense_421_525657:
dense_421_525659:,
batch_normalization_289_525662:,
batch_normalization_289_525664:,
batch_normalization_289_525666:,
batch_normalization_289_525668:"
dense_422_525681:
dense_422_525683:,
batch_normalization_290_525686:,
batch_normalization_290_525688:,
batch_normalization_290_525690:,
batch_normalization_290_525692:"
dense_423_525706:
dense_423_525708:
identity��/batch_normalization_288/StatefulPartitionedCall�/batch_normalization_289/StatefulPartitionedCall�/batch_normalization_290/StatefulPartitionedCall�!dense_420/StatefulPartitionedCall�!dense_421/StatefulPartitionedCall�!dense_422/StatefulPartitionedCall�!dense_423/StatefulPartitionedCall�
!dense_420/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_420_525633dense_420_525635*
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
E__inference_dense_420_layer_call_and_return_conditional_losses_525632�
/batch_normalization_288/StatefulPartitionedCallStatefulPartitionedCall*dense_420/StatefulPartitionedCall:output:0batch_normalization_288_525638batch_normalization_288_525640batch_normalization_288_525642batch_normalization_288_525644*
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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_525414�
!dense_421/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_288/StatefulPartitionedCall:output:0dense_421_525657dense_421_525659*
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
E__inference_dense_421_layer_call_and_return_conditional_losses_525656�
/batch_normalization_289/StatefulPartitionedCallStatefulPartitionedCall*dense_421/StatefulPartitionedCall:output:0batch_normalization_289_525662batch_normalization_289_525664batch_normalization_289_525666batch_normalization_289_525668*
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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_525494�
!dense_422/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_289/StatefulPartitionedCall:output:0dense_422_525681dense_422_525683*
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
E__inference_dense_422_layer_call_and_return_conditional_losses_525680�
/batch_normalization_290/StatefulPartitionedCallStatefulPartitionedCall*dense_422/StatefulPartitionedCall:output:0batch_normalization_290_525686batch_normalization_290_525688batch_normalization_290_525690batch_normalization_290_525692*
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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_525574�
!dense_423/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_290/StatefulPartitionedCall:output:0dense_423_525706dense_423_525708*
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
E__inference_dense_423_layer_call_and_return_conditional_losses_525705y
IdentityIdentity*dense_423/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp0^batch_normalization_288/StatefulPartitionedCall0^batch_normalization_289/StatefulPartitionedCall0^batch_normalization_290/StatefulPartitionedCall"^dense_420/StatefulPartitionedCall"^dense_421/StatefulPartitionedCall"^dense_422/StatefulPartitionedCall"^dense_423/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 2b
/batch_normalization_288/StatefulPartitionedCall/batch_normalization_288/StatefulPartitionedCall2b
/batch_normalization_289/StatefulPartitionedCall/batch_normalization_289/StatefulPartitionedCall2b
/batch_normalization_290/StatefulPartitionedCall/batch_normalization_290/StatefulPartitionedCall2F
!dense_420/StatefulPartitionedCall!dense_420/StatefulPartitionedCall2F
!dense_421/StatefulPartitionedCall!dense_421/StatefulPartitionedCall2F
!dense_422/StatefulPartitionedCall!dense_422/StatefulPartitionedCall2F
!dense_423/StatefulPartitionedCall!dense_423/StatefulPartitionedCall:&"
 
_user_specified_name525708:&"
 
_user_specified_name525706:&"
 
_user_specified_name525692:&"
 
_user_specified_name525690:&"
 
_user_specified_name525688:&"
 
_user_specified_name525686:&"
 
_user_specified_name525683:&"
 
_user_specified_name525681:&"
 
_user_specified_name525668:&"
 
_user_specified_name525666:&
"
 
_user_specified_name525664:&	"
 
_user_specified_name525662:&"
 
_user_specified_name525659:&"
 
_user_specified_name525657:&"
 
_user_specified_name525644:&"
 
_user_specified_name525642:&"
 
_user_specified_name525640:&"
 
_user_specified_name525638:&"
 
_user_specified_name525635:&"
 
_user_specified_name525633:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
0__inference_sparse_model_24_layer_call_fn_525808
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:

unknown_13:

unknown_14:

unknown_15:

unknown_16:

unknown_17:

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
GPU 2J 8� *T
fORM
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525712o
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
 
_user_specified_name525804:&"
 
_user_specified_name525802:&"
 
_user_specified_name525800:&"
 
_user_specified_name525798:&"
 
_user_specified_name525796:&"
 
_user_specified_name525794:&"
 
_user_specified_name525792:&"
 
_user_specified_name525790:&"
 
_user_specified_name525788:&"
 
_user_specified_name525786:&
"
 
_user_specified_name525784:&	"
 
_user_specified_name525782:&"
 
_user_specified_name525780:&"
 
_user_specified_name525778:&"
 
_user_specified_name525776:&"
 
_user_specified_name525774:&"
 
_user_specified_name525772:&"
 
_user_specified_name525770:&"
 
_user_specified_name525768:&"
 
_user_specified_name525766:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�
�
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526234

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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_525414

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
��
�
!__inference__wrapped_model_525380
input_1K
8sparse_model_24_dense_420_matmul_readvariableop_resource:	�G
9sparse_model_24_dense_420_biasadd_readvariableop_resource:R
Dsparse_model_24_batch_normalization_288_cast_readvariableop_resource:T
Fsparse_model_24_batch_normalization_288_cast_1_readvariableop_resource:T
Fsparse_model_24_batch_normalization_288_cast_2_readvariableop_resource:T
Fsparse_model_24_batch_normalization_288_cast_3_readvariableop_resource:J
8sparse_model_24_dense_421_matmul_readvariableop_resource:G
9sparse_model_24_dense_421_biasadd_readvariableop_resource:R
Dsparse_model_24_batch_normalization_289_cast_readvariableop_resource:T
Fsparse_model_24_batch_normalization_289_cast_1_readvariableop_resource:T
Fsparse_model_24_batch_normalization_289_cast_2_readvariableop_resource:T
Fsparse_model_24_batch_normalization_289_cast_3_readvariableop_resource:J
8sparse_model_24_dense_422_matmul_readvariableop_resource:G
9sparse_model_24_dense_422_biasadd_readvariableop_resource:R
Dsparse_model_24_batch_normalization_290_cast_readvariableop_resource:T
Fsparse_model_24_batch_normalization_290_cast_1_readvariableop_resource:T
Fsparse_model_24_batch_normalization_290_cast_2_readvariableop_resource:T
Fsparse_model_24_batch_normalization_290_cast_3_readvariableop_resource:J
8sparse_model_24_dense_423_matmul_readvariableop_resource:G
9sparse_model_24_dense_423_biasadd_readvariableop_resource:
identity��;sparse_model_24/batch_normalization_288/Cast/ReadVariableOp�=sparse_model_24/batch_normalization_288/Cast_1/ReadVariableOp�=sparse_model_24/batch_normalization_288/Cast_2/ReadVariableOp�=sparse_model_24/batch_normalization_288/Cast_3/ReadVariableOp�;sparse_model_24/batch_normalization_289/Cast/ReadVariableOp�=sparse_model_24/batch_normalization_289/Cast_1/ReadVariableOp�=sparse_model_24/batch_normalization_289/Cast_2/ReadVariableOp�=sparse_model_24/batch_normalization_289/Cast_3/ReadVariableOp�;sparse_model_24/batch_normalization_290/Cast/ReadVariableOp�=sparse_model_24/batch_normalization_290/Cast_1/ReadVariableOp�=sparse_model_24/batch_normalization_290/Cast_2/ReadVariableOp�=sparse_model_24/batch_normalization_290/Cast_3/ReadVariableOp�0sparse_model_24/dense_420/BiasAdd/ReadVariableOp�/sparse_model_24/dense_420/MatMul/ReadVariableOp�0sparse_model_24/dense_421/BiasAdd/ReadVariableOp�/sparse_model_24/dense_421/MatMul/ReadVariableOp�0sparse_model_24/dense_422/BiasAdd/ReadVariableOp�/sparse_model_24/dense_422/MatMul/ReadVariableOp�0sparse_model_24/dense_423/BiasAdd/ReadVariableOp�/sparse_model_24/dense_423/MatMul/ReadVariableOp�
/sparse_model_24/dense_420/MatMul/ReadVariableOpReadVariableOp8sparse_model_24_dense_420_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
 sparse_model_24/dense_420/MatMulMatMulinput_17sparse_model_24/dense_420/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
0sparse_model_24/dense_420/BiasAdd/ReadVariableOpReadVariableOp9sparse_model_24_dense_420_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
!sparse_model_24/dense_420/BiasAddBiasAdd*sparse_model_24/dense_420/MatMul:product:08sparse_model_24/dense_420/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
;sparse_model_24/batch_normalization_288/Cast/ReadVariableOpReadVariableOpDsparse_model_24_batch_normalization_288_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_288/Cast_1/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_288_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_288/Cast_2/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_288_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_288/Cast_3/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_288_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0|
7sparse_model_24/batch_normalization_288/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
5sparse_model_24/batch_normalization_288/batchnorm/addAddV2Esparse_model_24/batch_normalization_288/Cast_1/ReadVariableOp:value:0@sparse_model_24/batch_normalization_288/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_288/batchnorm/RsqrtRsqrt9sparse_model_24/batch_normalization_288/batchnorm/add:z:0*
T0*
_output_shapes
:�
5sparse_model_24/batch_normalization_288/batchnorm/mulMul;sparse_model_24/batch_normalization_288/batchnorm/Rsqrt:y:0Esparse_model_24/batch_normalization_288/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_288/batchnorm/mul_1Mul*sparse_model_24/dense_420/BiasAdd:output:09sparse_model_24/batch_normalization_288/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
7sparse_model_24/batch_normalization_288/batchnorm/mul_2MulCsparse_model_24/batch_normalization_288/Cast/ReadVariableOp:value:09sparse_model_24/batch_normalization_288/batchnorm/mul:z:0*
T0*
_output_shapes
:�
5sparse_model_24/batch_normalization_288/batchnorm/subSubEsparse_model_24/batch_normalization_288/Cast_2/ReadVariableOp:value:0;sparse_model_24/batch_normalization_288/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_288/batchnorm/add_1AddV2;sparse_model_24/batch_normalization_288/batchnorm/mul_1:z:09sparse_model_24/batch_normalization_288/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
/sparse_model_24/dense_421/MatMul/ReadVariableOpReadVariableOp8sparse_model_24_dense_421_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
 sparse_model_24/dense_421/MatMulMatMul;sparse_model_24/batch_normalization_288/batchnorm/add_1:z:07sparse_model_24/dense_421/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
0sparse_model_24/dense_421/BiasAdd/ReadVariableOpReadVariableOp9sparse_model_24_dense_421_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
!sparse_model_24/dense_421/BiasAddBiasAdd*sparse_model_24/dense_421/MatMul:product:08sparse_model_24/dense_421/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
;sparse_model_24/batch_normalization_289/Cast/ReadVariableOpReadVariableOpDsparse_model_24_batch_normalization_289_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_289/Cast_1/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_289_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_289/Cast_2/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_289_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_289/Cast_3/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_289_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0|
7sparse_model_24/batch_normalization_289/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
5sparse_model_24/batch_normalization_289/batchnorm/addAddV2Esparse_model_24/batch_normalization_289/Cast_1/ReadVariableOp:value:0@sparse_model_24/batch_normalization_289/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_289/batchnorm/RsqrtRsqrt9sparse_model_24/batch_normalization_289/batchnorm/add:z:0*
T0*
_output_shapes
:�
5sparse_model_24/batch_normalization_289/batchnorm/mulMul;sparse_model_24/batch_normalization_289/batchnorm/Rsqrt:y:0Esparse_model_24/batch_normalization_289/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_289/batchnorm/mul_1Mul*sparse_model_24/dense_421/BiasAdd:output:09sparse_model_24/batch_normalization_289/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
7sparse_model_24/batch_normalization_289/batchnorm/mul_2MulCsparse_model_24/batch_normalization_289/Cast/ReadVariableOp:value:09sparse_model_24/batch_normalization_289/batchnorm/mul:z:0*
T0*
_output_shapes
:�
5sparse_model_24/batch_normalization_289/batchnorm/subSubEsparse_model_24/batch_normalization_289/Cast_2/ReadVariableOp:value:0;sparse_model_24/batch_normalization_289/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_289/batchnorm/add_1AddV2;sparse_model_24/batch_normalization_289/batchnorm/mul_1:z:09sparse_model_24/batch_normalization_289/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
/sparse_model_24/dense_422/MatMul/ReadVariableOpReadVariableOp8sparse_model_24_dense_422_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
 sparse_model_24/dense_422/MatMulMatMul;sparse_model_24/batch_normalization_289/batchnorm/add_1:z:07sparse_model_24/dense_422/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
0sparse_model_24/dense_422/BiasAdd/ReadVariableOpReadVariableOp9sparse_model_24_dense_422_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
!sparse_model_24/dense_422/BiasAddBiasAdd*sparse_model_24/dense_422/MatMul:product:08sparse_model_24/dense_422/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
;sparse_model_24/batch_normalization_290/Cast/ReadVariableOpReadVariableOpDsparse_model_24_batch_normalization_290_cast_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_290/Cast_1/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_290_cast_1_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_290/Cast_2/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_290_cast_2_readvariableop_resource*
_output_shapes
:*
dtype0�
=sparse_model_24/batch_normalization_290/Cast_3/ReadVariableOpReadVariableOpFsparse_model_24_batch_normalization_290_cast_3_readvariableop_resource*
_output_shapes
:*
dtype0|
7sparse_model_24/batch_normalization_290/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o�:�
5sparse_model_24/batch_normalization_290/batchnorm/addAddV2Esparse_model_24/batch_normalization_290/Cast_1/ReadVariableOp:value:0@sparse_model_24/batch_normalization_290/batchnorm/add/y:output:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_290/batchnorm/RsqrtRsqrt9sparse_model_24/batch_normalization_290/batchnorm/add:z:0*
T0*
_output_shapes
:�
5sparse_model_24/batch_normalization_290/batchnorm/mulMul;sparse_model_24/batch_normalization_290/batchnorm/Rsqrt:y:0Esparse_model_24/batch_normalization_290/Cast_3/ReadVariableOp:value:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_290/batchnorm/mul_1Mul*sparse_model_24/dense_422/BiasAdd:output:09sparse_model_24/batch_normalization_290/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
7sparse_model_24/batch_normalization_290/batchnorm/mul_2MulCsparse_model_24/batch_normalization_290/Cast/ReadVariableOp:value:09sparse_model_24/batch_normalization_290/batchnorm/mul:z:0*
T0*
_output_shapes
:�
5sparse_model_24/batch_normalization_290/batchnorm/subSubEsparse_model_24/batch_normalization_290/Cast_2/ReadVariableOp:value:0;sparse_model_24/batch_normalization_290/batchnorm/mul_2:z:0*
T0*
_output_shapes
:�
7sparse_model_24/batch_normalization_290/batchnorm/add_1AddV2;sparse_model_24/batch_normalization_290/batchnorm/mul_1:z:09sparse_model_24/batch_normalization_290/batchnorm/sub:z:0*
T0*'
_output_shapes
:����������
/sparse_model_24/dense_423/MatMul/ReadVariableOpReadVariableOp8sparse_model_24_dense_423_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
 sparse_model_24/dense_423/MatMulMatMul;sparse_model_24/batch_normalization_290/batchnorm/add_1:z:07sparse_model_24/dense_423/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
0sparse_model_24/dense_423/BiasAdd/ReadVariableOpReadVariableOp9sparse_model_24_dense_423_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
!sparse_model_24/dense_423/BiasAddBiasAdd*sparse_model_24/dense_423/MatMul:product:08sparse_model_24/dense_423/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
!sparse_model_24/dense_423/SigmoidSigmoid*sparse_model_24/dense_423/BiasAdd:output:0*
T0*'
_output_shapes
:���������t
IdentityIdentity%sparse_model_24/dense_423/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:����������	
NoOpNoOp<^sparse_model_24/batch_normalization_288/Cast/ReadVariableOp>^sparse_model_24/batch_normalization_288/Cast_1/ReadVariableOp>^sparse_model_24/batch_normalization_288/Cast_2/ReadVariableOp>^sparse_model_24/batch_normalization_288/Cast_3/ReadVariableOp<^sparse_model_24/batch_normalization_289/Cast/ReadVariableOp>^sparse_model_24/batch_normalization_289/Cast_1/ReadVariableOp>^sparse_model_24/batch_normalization_289/Cast_2/ReadVariableOp>^sparse_model_24/batch_normalization_289/Cast_3/ReadVariableOp<^sparse_model_24/batch_normalization_290/Cast/ReadVariableOp>^sparse_model_24/batch_normalization_290/Cast_1/ReadVariableOp>^sparse_model_24/batch_normalization_290/Cast_2/ReadVariableOp>^sparse_model_24/batch_normalization_290/Cast_3/ReadVariableOp1^sparse_model_24/dense_420/BiasAdd/ReadVariableOp0^sparse_model_24/dense_420/MatMul/ReadVariableOp1^sparse_model_24/dense_421/BiasAdd/ReadVariableOp0^sparse_model_24/dense_421/MatMul/ReadVariableOp1^sparse_model_24/dense_422/BiasAdd/ReadVariableOp0^sparse_model_24/dense_422/MatMul/ReadVariableOp1^sparse_model_24/dense_423/BiasAdd/ReadVariableOp0^sparse_model_24/dense_423/MatMul/ReadVariableOp*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 2z
;sparse_model_24/batch_normalization_288/Cast/ReadVariableOp;sparse_model_24/batch_normalization_288/Cast/ReadVariableOp2~
=sparse_model_24/batch_normalization_288/Cast_1/ReadVariableOp=sparse_model_24/batch_normalization_288/Cast_1/ReadVariableOp2~
=sparse_model_24/batch_normalization_288/Cast_2/ReadVariableOp=sparse_model_24/batch_normalization_288/Cast_2/ReadVariableOp2~
=sparse_model_24/batch_normalization_288/Cast_3/ReadVariableOp=sparse_model_24/batch_normalization_288/Cast_3/ReadVariableOp2z
;sparse_model_24/batch_normalization_289/Cast/ReadVariableOp;sparse_model_24/batch_normalization_289/Cast/ReadVariableOp2~
=sparse_model_24/batch_normalization_289/Cast_1/ReadVariableOp=sparse_model_24/batch_normalization_289/Cast_1/ReadVariableOp2~
=sparse_model_24/batch_normalization_289/Cast_2/ReadVariableOp=sparse_model_24/batch_normalization_289/Cast_2/ReadVariableOp2~
=sparse_model_24/batch_normalization_289/Cast_3/ReadVariableOp=sparse_model_24/batch_normalization_289/Cast_3/ReadVariableOp2z
;sparse_model_24/batch_normalization_290/Cast/ReadVariableOp;sparse_model_24/batch_normalization_290/Cast/ReadVariableOp2~
=sparse_model_24/batch_normalization_290/Cast_1/ReadVariableOp=sparse_model_24/batch_normalization_290/Cast_1/ReadVariableOp2~
=sparse_model_24/batch_normalization_290/Cast_2/ReadVariableOp=sparse_model_24/batch_normalization_290/Cast_2/ReadVariableOp2~
=sparse_model_24/batch_normalization_290/Cast_3/ReadVariableOp=sparse_model_24/batch_normalization_290/Cast_3/ReadVariableOp2d
0sparse_model_24/dense_420/BiasAdd/ReadVariableOp0sparse_model_24/dense_420/BiasAdd/ReadVariableOp2b
/sparse_model_24/dense_420/MatMul/ReadVariableOp/sparse_model_24/dense_420/MatMul/ReadVariableOp2d
0sparse_model_24/dense_421/BiasAdd/ReadVariableOp0sparse_model_24/dense_421/BiasAdd/ReadVariableOp2b
/sparse_model_24/dense_421/MatMul/ReadVariableOp/sparse_model_24/dense_421/MatMul/ReadVariableOp2d
0sparse_model_24/dense_422/BiasAdd/ReadVariableOp0sparse_model_24/dense_422/BiasAdd/ReadVariableOp2b
/sparse_model_24/dense_422/MatMul/ReadVariableOp/sparse_model_24/dense_422/MatMul/ReadVariableOp2d
0sparse_model_24/dense_423/BiasAdd/ReadVariableOp0sparse_model_24/dense_423/BiasAdd/ReadVariableOp2b
/sparse_model_24/dense_423/MatMul/ReadVariableOp/sparse_model_24/dense_423/MatMul/ReadVariableOp:($
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
8__inference_batch_normalization_290_layer_call_fn_526180

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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_525594o
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
 
_user_specified_name526176:&"
 
_user_specified_name526174:&"
 
_user_specified_name526172:&"
 
_user_specified_name526170:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
*__inference_dense_421_layer_call_fn_526045

inputs
unknown:
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
E__inference_dense_421_layer_call_and_return_conditional_losses_525656o
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
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:&"
 
_user_specified_name526041:&"
 
_user_specified_name526039:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526135

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
�
�
$__inference_signature_wrapper_525937
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:

unknown_13:

unknown_14:

unknown_15:

unknown_16:

unknown_17:

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
!__inference__wrapped_model_525380o
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
 
_user_specified_name525933:&"
 
_user_specified_name525931:&"
 
_user_specified_name525929:&"
 
_user_specified_name525927:&"
 
_user_specified_name525925:&"
 
_user_specified_name525923:&"
 
_user_specified_name525921:&"
 
_user_specified_name525919:&"
 
_user_specified_name525917:&"
 
_user_specified_name525915:&
"
 
_user_specified_name525913:&	"
 
_user_specified_name525911:&"
 
_user_specified_name525909:&"
 
_user_specified_name525907:&"
 
_user_specified_name525905:&"
 
_user_specified_name525903:&"
 
_user_specified_name525901:&"
 
_user_specified_name525899:&"
 
_user_specified_name525897:&"
 
_user_specified_name525895:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�%
�
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_525574

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
8__inference_batch_normalization_289_layer_call_fn_526068

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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_525494o
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
 
_user_specified_name526064:&"
 
_user_specified_name526062:&"
 
_user_specified_name526060:&"
 
_user_specified_name526058:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�%
�
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_525494

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
E__inference_dense_421_layer_call_and_return_conditional_losses_525656

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_525434

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
E__inference_dense_420_layer_call_and_return_conditional_losses_525632

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
E__inference_dense_420_layer_call_and_return_conditional_losses_525956

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
�
�
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_525594

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
�)
"__inference__traced_restore_526771
file_prefixD
1assignvariableop_sparse_model_24_dense_420_kernel:	�?
1assignvariableop_1_sparse_model_24_dense_420_bias:N
@assignvariableop_2_sparse_model_24_batch_normalization_288_gamma:M
?assignvariableop_3_sparse_model_24_batch_normalization_288_beta:T
Fassignvariableop_4_sparse_model_24_batch_normalization_288_moving_mean:X
Jassignvariableop_5_sparse_model_24_batch_normalization_288_moving_variance:E
3assignvariableop_6_sparse_model_24_dense_421_kernel:?
1assignvariableop_7_sparse_model_24_dense_421_bias:N
@assignvariableop_8_sparse_model_24_batch_normalization_289_gamma:M
?assignvariableop_9_sparse_model_24_batch_normalization_289_beta:U
Gassignvariableop_10_sparse_model_24_batch_normalization_289_moving_mean:Y
Kassignvariableop_11_sparse_model_24_batch_normalization_289_moving_variance:F
4assignvariableop_12_sparse_model_24_dense_422_kernel:@
2assignvariableop_13_sparse_model_24_dense_422_bias:O
Aassignvariableop_14_sparse_model_24_batch_normalization_290_gamma:N
@assignvariableop_15_sparse_model_24_batch_normalization_290_beta:U
Gassignvariableop_16_sparse_model_24_batch_normalization_290_moving_mean:Y
Kassignvariableop_17_sparse_model_24_batch_normalization_290_moving_variance:F
4assignvariableop_18_sparse_model_24_dense_423_kernel:@
2assignvariableop_19_sparse_model_24_dense_423_bias:'
assignvariableop_20_iteration:	 +
!assignvariableop_21_learning_rate: N
;assignvariableop_22_adam_m_sparse_model_24_dense_420_kernel:	�N
;assignvariableop_23_adam_v_sparse_model_24_dense_420_kernel:	�G
9assignvariableop_24_adam_m_sparse_model_24_dense_420_bias:G
9assignvariableop_25_adam_v_sparse_model_24_dense_420_bias:V
Hassignvariableop_26_adam_m_sparse_model_24_batch_normalization_288_gamma:V
Hassignvariableop_27_adam_v_sparse_model_24_batch_normalization_288_gamma:U
Gassignvariableop_28_adam_m_sparse_model_24_batch_normalization_288_beta:U
Gassignvariableop_29_adam_v_sparse_model_24_batch_normalization_288_beta:M
;assignvariableop_30_adam_m_sparse_model_24_dense_421_kernel:M
;assignvariableop_31_adam_v_sparse_model_24_dense_421_kernel:G
9assignvariableop_32_adam_m_sparse_model_24_dense_421_bias:G
9assignvariableop_33_adam_v_sparse_model_24_dense_421_bias:V
Hassignvariableop_34_adam_m_sparse_model_24_batch_normalization_289_gamma:V
Hassignvariableop_35_adam_v_sparse_model_24_batch_normalization_289_gamma:U
Gassignvariableop_36_adam_m_sparse_model_24_batch_normalization_289_beta:U
Gassignvariableop_37_adam_v_sparse_model_24_batch_normalization_289_beta:M
;assignvariableop_38_adam_m_sparse_model_24_dense_422_kernel:M
;assignvariableop_39_adam_v_sparse_model_24_dense_422_kernel:G
9assignvariableop_40_adam_m_sparse_model_24_dense_422_bias:G
9assignvariableop_41_adam_v_sparse_model_24_dense_422_bias:V
Hassignvariableop_42_adam_m_sparse_model_24_batch_normalization_290_gamma:V
Hassignvariableop_43_adam_v_sparse_model_24_batch_normalization_290_gamma:U
Gassignvariableop_44_adam_m_sparse_model_24_batch_normalization_290_beta:U
Gassignvariableop_45_adam_v_sparse_model_24_batch_normalization_290_beta:M
;assignvariableop_46_adam_m_sparse_model_24_dense_423_kernel:M
;assignvariableop_47_adam_v_sparse_model_24_dense_423_kernel:G
9assignvariableop_48_adam_m_sparse_model_24_dense_423_bias:G
9assignvariableop_49_adam_v_sparse_model_24_dense_423_bias:%
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
AssignVariableOpAssignVariableOp1assignvariableop_sparse_model_24_dense_420_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp1assignvariableop_1_sparse_model_24_dense_420_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOp@assignvariableop_2_sparse_model_24_batch_normalization_288_gammaIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOp?assignvariableop_3_sparse_model_24_batch_normalization_288_betaIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpFassignvariableop_4_sparse_model_24_batch_normalization_288_moving_meanIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpJassignvariableop_5_sparse_model_24_batch_normalization_288_moving_varianceIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp3assignvariableop_6_sparse_model_24_dense_421_kernelIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp1assignvariableop_7_sparse_model_24_dense_421_biasIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOp@assignvariableop_8_sparse_model_24_batch_normalization_289_gammaIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp?assignvariableop_9_sparse_model_24_batch_normalization_289_betaIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOpGassignvariableop_10_sparse_model_24_batch_normalization_289_moving_meanIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpKassignvariableop_11_sparse_model_24_batch_normalization_289_moving_varianceIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp4assignvariableop_12_sparse_model_24_dense_422_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp2assignvariableop_13_sparse_model_24_dense_422_biasIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOpAassignvariableop_14_sparse_model_24_batch_normalization_290_gammaIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOp@assignvariableop_15_sparse_model_24_batch_normalization_290_betaIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOpGassignvariableop_16_sparse_model_24_batch_normalization_290_moving_meanIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOpKassignvariableop_17_sparse_model_24_batch_normalization_290_moving_varianceIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp4assignvariableop_18_sparse_model_24_dense_423_kernelIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp2assignvariableop_19_sparse_model_24_dense_423_biasIdentity_19:output:0"/device:CPU:0*&
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
AssignVariableOp_22AssignVariableOp;assignvariableop_22_adam_m_sparse_model_24_dense_420_kernelIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOp;assignvariableop_23_adam_v_sparse_model_24_dense_420_kernelIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOp9assignvariableop_24_adam_m_sparse_model_24_dense_420_biasIdentity_24:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOp9assignvariableop_25_adam_v_sparse_model_24_dense_420_biasIdentity_25:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpHassignvariableop_26_adam_m_sparse_model_24_batch_normalization_288_gammaIdentity_26:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOpHassignvariableop_27_adam_v_sparse_model_24_batch_normalization_288_gammaIdentity_27:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpGassignvariableop_28_adam_m_sparse_model_24_batch_normalization_288_betaIdentity_28:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpGassignvariableop_29_adam_v_sparse_model_24_batch_normalization_288_betaIdentity_29:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOp;assignvariableop_30_adam_m_sparse_model_24_dense_421_kernelIdentity_30:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOp;assignvariableop_31_adam_v_sparse_model_24_dense_421_kernelIdentity_31:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOp9assignvariableop_32_adam_m_sparse_model_24_dense_421_biasIdentity_32:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOp9assignvariableop_33_adam_v_sparse_model_24_dense_421_biasIdentity_33:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOpHassignvariableop_34_adam_m_sparse_model_24_batch_normalization_289_gammaIdentity_34:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOpHassignvariableop_35_adam_v_sparse_model_24_batch_normalization_289_gammaIdentity_35:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_36AssignVariableOpGassignvariableop_36_adam_m_sparse_model_24_batch_normalization_289_betaIdentity_36:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_37AssignVariableOpGassignvariableop_37_adam_v_sparse_model_24_batch_normalization_289_betaIdentity_37:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_38AssignVariableOp;assignvariableop_38_adam_m_sparse_model_24_dense_422_kernelIdentity_38:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_39AssignVariableOp;assignvariableop_39_adam_v_sparse_model_24_dense_422_kernelIdentity_39:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_40AssignVariableOp9assignvariableop_40_adam_m_sparse_model_24_dense_422_biasIdentity_40:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_41AssignVariableOp9assignvariableop_41_adam_v_sparse_model_24_dense_422_biasIdentity_41:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_42AssignVariableOpHassignvariableop_42_adam_m_sparse_model_24_batch_normalization_290_gammaIdentity_42:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_43AssignVariableOpHassignvariableop_43_adam_v_sparse_model_24_batch_normalization_290_gammaIdentity_43:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_44AssignVariableOpGassignvariableop_44_adam_m_sparse_model_24_batch_normalization_290_betaIdentity_44:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_45AssignVariableOpGassignvariableop_45_adam_v_sparse_model_24_batch_normalization_290_betaIdentity_45:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_46AssignVariableOp;assignvariableop_46_adam_m_sparse_model_24_dense_423_kernelIdentity_46:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_47AssignVariableOp;assignvariableop_47_adam_v_sparse_model_24_dense_423_kernelIdentity_47:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_48AssignVariableOp9assignvariableop_48_adam_m_sparse_model_24_dense_423_biasIdentity_48:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_49AssignVariableOp9assignvariableop_49_adam_v_sparse_model_24_dense_423_biasIdentity_49:output:0"/device:CPU:0*&
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
_user_specified_name	total_1:E2A
?
_user_specified_name'%Adam/v/sparse_model_24/dense_423/bias:E1A
?
_user_specified_name'%Adam/m/sparse_model_24/dense_423/bias:G0C
A
_user_specified_name)'Adam/v/sparse_model_24/dense_423/kernel:G/C
A
_user_specified_name)'Adam/m/sparse_model_24/dense_423/kernel:S.O
M
_user_specified_name53Adam/v/sparse_model_24/batch_normalization_290/beta:S-O
M
_user_specified_name53Adam/m/sparse_model_24/batch_normalization_290/beta:T,P
N
_user_specified_name64Adam/v/sparse_model_24/batch_normalization_290/gamma:T+P
N
_user_specified_name64Adam/m/sparse_model_24/batch_normalization_290/gamma:E*A
?
_user_specified_name'%Adam/v/sparse_model_24/dense_422/bias:E)A
?
_user_specified_name'%Adam/m/sparse_model_24/dense_422/bias:G(C
A
_user_specified_name)'Adam/v/sparse_model_24/dense_422/kernel:G'C
A
_user_specified_name)'Adam/m/sparse_model_24/dense_422/kernel:S&O
M
_user_specified_name53Adam/v/sparse_model_24/batch_normalization_289/beta:S%O
M
_user_specified_name53Adam/m/sparse_model_24/batch_normalization_289/beta:T$P
N
_user_specified_name64Adam/v/sparse_model_24/batch_normalization_289/gamma:T#P
N
_user_specified_name64Adam/m/sparse_model_24/batch_normalization_289/gamma:E"A
?
_user_specified_name'%Adam/v/sparse_model_24/dense_421/bias:E!A
?
_user_specified_name'%Adam/m/sparse_model_24/dense_421/bias:G C
A
_user_specified_name)'Adam/v/sparse_model_24/dense_421/kernel:GC
A
_user_specified_name)'Adam/m/sparse_model_24/dense_421/kernel:SO
M
_user_specified_name53Adam/v/sparse_model_24/batch_normalization_288/beta:SO
M
_user_specified_name53Adam/m/sparse_model_24/batch_normalization_288/beta:TP
N
_user_specified_name64Adam/v/sparse_model_24/batch_normalization_288/gamma:TP
N
_user_specified_name64Adam/m/sparse_model_24/batch_normalization_288/gamma:EA
?
_user_specified_name'%Adam/v/sparse_model_24/dense_420/bias:EA
?
_user_specified_name'%Adam/m/sparse_model_24/dense_420/bias:GC
A
_user_specified_name)'Adam/v/sparse_model_24/dense_420/kernel:GC
A
_user_specified_name)'Adam/m/sparse_model_24/dense_420/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:>:
8
_user_specified_name sparse_model_24/dense_423/bias:@<
:
_user_specified_name" sparse_model_24/dense_423/kernel:WS
Q
_user_specified_name97sparse_model_24/batch_normalization_290/moving_variance:SO
M
_user_specified_name53sparse_model_24/batch_normalization_290/moving_mean:LH
F
_user_specified_name.,sparse_model_24/batch_normalization_290/beta:MI
G
_user_specified_name/-sparse_model_24/batch_normalization_290/gamma:>:
8
_user_specified_name sparse_model_24/dense_422/bias:@<
:
_user_specified_name" sparse_model_24/dense_422/kernel:WS
Q
_user_specified_name97sparse_model_24/batch_normalization_289/moving_variance:SO
M
_user_specified_name53sparse_model_24/batch_normalization_289/moving_mean:L
H
F
_user_specified_name.,sparse_model_24/batch_normalization_289/beta:M	I
G
_user_specified_name/-sparse_model_24/batch_normalization_289/gamma:>:
8
_user_specified_name sparse_model_24/dense_421/bias:@<
:
_user_specified_name" sparse_model_24/dense_421/kernel:WS
Q
_user_specified_name97sparse_model_24/batch_normalization_288/moving_variance:SO
M
_user_specified_name53sparse_model_24/batch_normalization_288/moving_mean:LH
F
_user_specified_name.,sparse_model_24/batch_normalization_288/beta:MI
G
_user_specified_name/-sparse_model_24/batch_normalization_288/gamma:>:
8
_user_specified_name sparse_model_24/dense_420/bias:@<
:
_user_specified_name" sparse_model_24/dense_420/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�3
�	
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525763
input_1#
dense_420_525715:	�
dense_420_525717:,
batch_normalization_288_525720:,
batch_normalization_288_525722:,
batch_normalization_288_525724:,
batch_normalization_288_525726:"
dense_421_525729:
dense_421_525731:,
batch_normalization_289_525734:,
batch_normalization_289_525736:,
batch_normalization_289_525738:,
batch_normalization_289_525740:"
dense_422_525743:
dense_422_525745:,
batch_normalization_290_525748:,
batch_normalization_290_525750:,
batch_normalization_290_525752:,
batch_normalization_290_525754:"
dense_423_525757:
dense_423_525759:
identity��/batch_normalization_288/StatefulPartitionedCall�/batch_normalization_289/StatefulPartitionedCall�/batch_normalization_290/StatefulPartitionedCall�!dense_420/StatefulPartitionedCall�!dense_421/StatefulPartitionedCall�!dense_422/StatefulPartitionedCall�!dense_423/StatefulPartitionedCall�
!dense_420/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_420_525715dense_420_525717*
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
E__inference_dense_420_layer_call_and_return_conditional_losses_525632�
/batch_normalization_288/StatefulPartitionedCallStatefulPartitionedCall*dense_420/StatefulPartitionedCall:output:0batch_normalization_288_525720batch_normalization_288_525722batch_normalization_288_525724batch_normalization_288_525726*
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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_525434�
!dense_421/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_288/StatefulPartitionedCall:output:0dense_421_525729dense_421_525731*
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
E__inference_dense_421_layer_call_and_return_conditional_losses_525656�
/batch_normalization_289/StatefulPartitionedCallStatefulPartitionedCall*dense_421/StatefulPartitionedCall:output:0batch_normalization_289_525734batch_normalization_289_525736batch_normalization_289_525738batch_normalization_289_525740*
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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_525514�
!dense_422/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_289/StatefulPartitionedCall:output:0dense_422_525743dense_422_525745*
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
E__inference_dense_422_layer_call_and_return_conditional_losses_525680�
/batch_normalization_290/StatefulPartitionedCallStatefulPartitionedCall*dense_422/StatefulPartitionedCall:output:0batch_normalization_290_525748batch_normalization_290_525750batch_normalization_290_525752batch_normalization_290_525754*
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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_525594�
!dense_423/StatefulPartitionedCallStatefulPartitionedCall8batch_normalization_290/StatefulPartitionedCall:output:0dense_423_525757dense_423_525759*
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
E__inference_dense_423_layer_call_and_return_conditional_losses_525705y
IdentityIdentity*dense_423/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp0^batch_normalization_288/StatefulPartitionedCall0^batch_normalization_289/StatefulPartitionedCall0^batch_normalization_290/StatefulPartitionedCall"^dense_420/StatefulPartitionedCall"^dense_421/StatefulPartitionedCall"^dense_422/StatefulPartitionedCall"^dense_423/StatefulPartitionedCall*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*O
_input_shapes>
<:����������: : : : : : : : : : : : : : : : : : : : 2b
/batch_normalization_288/StatefulPartitionedCall/batch_normalization_288/StatefulPartitionedCall2b
/batch_normalization_289/StatefulPartitionedCall/batch_normalization_289/StatefulPartitionedCall2b
/batch_normalization_290/StatefulPartitionedCall/batch_normalization_290/StatefulPartitionedCall2F
!dense_420/StatefulPartitionedCall!dense_420/StatefulPartitionedCall2F
!dense_421/StatefulPartitionedCall!dense_421/StatefulPartitionedCall2F
!dense_422/StatefulPartitionedCall!dense_422/StatefulPartitionedCall2F
!dense_423/StatefulPartitionedCall!dense_423/StatefulPartitionedCall:&"
 
_user_specified_name525759:&"
 
_user_specified_name525757:&"
 
_user_specified_name525754:&"
 
_user_specified_name525752:&"
 
_user_specified_name525750:&"
 
_user_specified_name525748:&"
 
_user_specified_name525745:&"
 
_user_specified_name525743:&"
 
_user_specified_name525740:&"
 
_user_specified_name525738:&
"
 
_user_specified_name525736:&	"
 
_user_specified_name525734:&"
 
_user_specified_name525731:&"
 
_user_specified_name525729:&"
 
_user_specified_name525726:&"
 
_user_specified_name525724:&"
 
_user_specified_name525722:&"
 
_user_specified_name525720:&"
 
_user_specified_name525717:&"
 
_user_specified_name525715:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
�

�
E__inference_dense_423_layer_call_and_return_conditional_losses_525705

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
�	
�
E__inference_dense_422_layer_call_and_return_conditional_losses_526154

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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_525514

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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526214

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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526036

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
8__inference_batch_normalization_289_layer_call_fn_526081

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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_525514o
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
 
_user_specified_name526077:&"
 
_user_specified_name526075:&"
 
_user_specified_name526073:&"
 
_user_specified_name526071:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_421_layer_call_and_return_conditional_losses_526055

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
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
�

�
E__inference_dense_423_layer_call_and_return_conditional_losses_526254

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
�	
�
8__inference_batch_normalization_288_layer_call_fn_525982

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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_525434o
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
 
_user_specified_name525978:&"
 
_user_specified_name525976:&"
 
_user_specified_name525974:&"
 
_user_specified_name525972:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
8__inference_batch_normalization_288_layer_call_fn_525969

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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_525414o
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
 
_user_specified_name525965:&"
 
_user_specified_name525963:&"
 
_user_specified_name525961:&"
 
_user_specified_name525959:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
8__inference_batch_normalization_290_layer_call_fn_526167

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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_525574o
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
 
_user_specified_name526163:&"
 
_user_specified_name526161:&"
 
_user_specified_name526159:&"
 
_user_specified_name526157:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
*__inference_dense_420_layer_call_fn_525946

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
E__inference_dense_420_layer_call_and_return_conditional_losses_525632o
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
 
_user_specified_name525942:&"
 
_user_specified_name525940:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
*__inference_dense_422_layer_call_fn_526144

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
E__inference_dense_422_layer_call_and_return_conditional_losses_525680o
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
 
_user_specified_name526140:&"
 
_user_specified_name526138:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�%
�
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526016

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
�
�
*__inference_dense_423_layer_call_fn_526243

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
E__inference_dense_423_layer_call_and_return_conditional_losses_525705o
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
 
_user_specified_name526239:&"
 
_user_specified_name526237:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
E__inference_dense_422_layer_call_and_return_conditional_losses_525680

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
�
0__inference_sparse_model_24_layer_call_fn_525853
input_1
unknown:	�
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3:
	unknown_4:
	unknown_5:
	unknown_6:
	unknown_7:
	unknown_8:
	unknown_9:

unknown_10:

unknown_11:

unknown_12:

unknown_13:

unknown_14:

unknown_15:

unknown_16:

unknown_17:

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
GPU 2J 8� *T
fORM
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525763o
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
 
_user_specified_name525849:&"
 
_user_specified_name525847:&"
 
_user_specified_name525845:&"
 
_user_specified_name525843:&"
 
_user_specified_name525841:&"
 
_user_specified_name525839:&"
 
_user_specified_name525837:&"
 
_user_specified_name525835:&"
 
_user_specified_name525833:&"
 
_user_specified_name525831:&
"
 
_user_specified_name525829:&	"
 
_user_specified_name525827:&"
 
_user_specified_name525825:&"
 
_user_specified_name525823:&"
 
_user_specified_name525821:&"
 
_user_specified_name525819:&"
 
_user_specified_name525817:&"
 
_user_specified_name525815:&"
 
_user_specified_name525813:&"
 
_user_specified_name525811:Q M
(
_output_shapes
:����������
!
_user_specified_name	input_1
˵
�9
__inference__traced_save_526600
file_prefixJ
7read_disablecopyonread_sparse_model_24_dense_420_kernel:	�E
7read_1_disablecopyonread_sparse_model_24_dense_420_bias:T
Fread_2_disablecopyonread_sparse_model_24_batch_normalization_288_gamma:S
Eread_3_disablecopyonread_sparse_model_24_batch_normalization_288_beta:Z
Lread_4_disablecopyonread_sparse_model_24_batch_normalization_288_moving_mean:^
Pread_5_disablecopyonread_sparse_model_24_batch_normalization_288_moving_variance:K
9read_6_disablecopyonread_sparse_model_24_dense_421_kernel:E
7read_7_disablecopyonread_sparse_model_24_dense_421_bias:T
Fread_8_disablecopyonread_sparse_model_24_batch_normalization_289_gamma:S
Eread_9_disablecopyonread_sparse_model_24_batch_normalization_289_beta:[
Mread_10_disablecopyonread_sparse_model_24_batch_normalization_289_moving_mean:_
Qread_11_disablecopyonread_sparse_model_24_batch_normalization_289_moving_variance:L
:read_12_disablecopyonread_sparse_model_24_dense_422_kernel:F
8read_13_disablecopyonread_sparse_model_24_dense_422_bias:U
Gread_14_disablecopyonread_sparse_model_24_batch_normalization_290_gamma:T
Fread_15_disablecopyonread_sparse_model_24_batch_normalization_290_beta:[
Mread_16_disablecopyonread_sparse_model_24_batch_normalization_290_moving_mean:_
Qread_17_disablecopyonread_sparse_model_24_batch_normalization_290_moving_variance:L
:read_18_disablecopyonread_sparse_model_24_dense_423_kernel:F
8read_19_disablecopyonread_sparse_model_24_dense_423_bias:-
#read_20_disablecopyonread_iteration:	 1
'read_21_disablecopyonread_learning_rate: T
Aread_22_disablecopyonread_adam_m_sparse_model_24_dense_420_kernel:	�T
Aread_23_disablecopyonread_adam_v_sparse_model_24_dense_420_kernel:	�M
?read_24_disablecopyonread_adam_m_sparse_model_24_dense_420_bias:M
?read_25_disablecopyonread_adam_v_sparse_model_24_dense_420_bias:\
Nread_26_disablecopyonread_adam_m_sparse_model_24_batch_normalization_288_gamma:\
Nread_27_disablecopyonread_adam_v_sparse_model_24_batch_normalization_288_gamma:[
Mread_28_disablecopyonread_adam_m_sparse_model_24_batch_normalization_288_beta:[
Mread_29_disablecopyonread_adam_v_sparse_model_24_batch_normalization_288_beta:S
Aread_30_disablecopyonread_adam_m_sparse_model_24_dense_421_kernel:S
Aread_31_disablecopyonread_adam_v_sparse_model_24_dense_421_kernel:M
?read_32_disablecopyonread_adam_m_sparse_model_24_dense_421_bias:M
?read_33_disablecopyonread_adam_v_sparse_model_24_dense_421_bias:\
Nread_34_disablecopyonread_adam_m_sparse_model_24_batch_normalization_289_gamma:\
Nread_35_disablecopyonread_adam_v_sparse_model_24_batch_normalization_289_gamma:[
Mread_36_disablecopyonread_adam_m_sparse_model_24_batch_normalization_289_beta:[
Mread_37_disablecopyonread_adam_v_sparse_model_24_batch_normalization_289_beta:S
Aread_38_disablecopyonread_adam_m_sparse_model_24_dense_422_kernel:S
Aread_39_disablecopyonread_adam_v_sparse_model_24_dense_422_kernel:M
?read_40_disablecopyonread_adam_m_sparse_model_24_dense_422_bias:M
?read_41_disablecopyonread_adam_v_sparse_model_24_dense_422_bias:\
Nread_42_disablecopyonread_adam_m_sparse_model_24_batch_normalization_290_gamma:\
Nread_43_disablecopyonread_adam_v_sparse_model_24_batch_normalization_290_gamma:[
Mread_44_disablecopyonread_adam_m_sparse_model_24_batch_normalization_290_beta:[
Mread_45_disablecopyonread_adam_v_sparse_model_24_batch_normalization_290_beta:S
Aread_46_disablecopyonread_adam_m_sparse_model_24_dense_423_kernel:S
Aread_47_disablecopyonread_adam_v_sparse_model_24_dense_423_kernel:M
?read_48_disablecopyonread_adam_m_sparse_model_24_dense_423_bias:M
?read_49_disablecopyonread_adam_v_sparse_model_24_dense_423_bias:+
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
Read/DisableCopyOnReadDisableCopyOnRead7read_disablecopyonread_sparse_model_24_dense_420_kernel"/device:CPU:0*
_output_shapes
 �
Read/ReadVariableOpReadVariableOp7read_disablecopyonread_sparse_model_24_dense_420_kernel^Read/DisableCopyOnRead"/device:CPU:0*
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
Read_1/DisableCopyOnReadDisableCopyOnRead7read_1_disablecopyonread_sparse_model_24_dense_420_bias"/device:CPU:0*
_output_shapes
 �
Read_1/ReadVariableOpReadVariableOp7read_1_disablecopyonread_sparse_model_24_dense_420_bias^Read_1/DisableCopyOnRead"/device:CPU:0*
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
Read_2/DisableCopyOnReadDisableCopyOnReadFread_2_disablecopyonread_sparse_model_24_batch_normalization_288_gamma"/device:CPU:0*
_output_shapes
 �
Read_2/ReadVariableOpReadVariableOpFread_2_disablecopyonread_sparse_model_24_batch_normalization_288_gamma^Read_2/DisableCopyOnRead"/device:CPU:0*
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
Read_3/DisableCopyOnReadDisableCopyOnReadEread_3_disablecopyonread_sparse_model_24_batch_normalization_288_beta"/device:CPU:0*
_output_shapes
 �
Read_3/ReadVariableOpReadVariableOpEread_3_disablecopyonread_sparse_model_24_batch_normalization_288_beta^Read_3/DisableCopyOnRead"/device:CPU:0*
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
Read_4/DisableCopyOnReadDisableCopyOnReadLread_4_disablecopyonread_sparse_model_24_batch_normalization_288_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_4/ReadVariableOpReadVariableOpLread_4_disablecopyonread_sparse_model_24_batch_normalization_288_moving_mean^Read_4/DisableCopyOnRead"/device:CPU:0*
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
Read_5/DisableCopyOnReadDisableCopyOnReadPread_5_disablecopyonread_sparse_model_24_batch_normalization_288_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_5/ReadVariableOpReadVariableOpPread_5_disablecopyonread_sparse_model_24_batch_normalization_288_moving_variance^Read_5/DisableCopyOnRead"/device:CPU:0*
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
Read_6/DisableCopyOnReadDisableCopyOnRead9read_6_disablecopyonread_sparse_model_24_dense_421_kernel"/device:CPU:0*
_output_shapes
 �
Read_6/ReadVariableOpReadVariableOp9read_6_disablecopyonread_sparse_model_24_dense_421_kernel^Read_6/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0n
Identity_12IdentityRead_6/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_13IdentityIdentity_12:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_7/DisableCopyOnReadDisableCopyOnRead7read_7_disablecopyonread_sparse_model_24_dense_421_bias"/device:CPU:0*
_output_shapes
 �
Read_7/ReadVariableOpReadVariableOp7read_7_disablecopyonread_sparse_model_24_dense_421_bias^Read_7/DisableCopyOnRead"/device:CPU:0*
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
Read_8/DisableCopyOnReadDisableCopyOnReadFread_8_disablecopyonread_sparse_model_24_batch_normalization_289_gamma"/device:CPU:0*
_output_shapes
 �
Read_8/ReadVariableOpReadVariableOpFread_8_disablecopyonread_sparse_model_24_batch_normalization_289_gamma^Read_8/DisableCopyOnRead"/device:CPU:0*
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
Read_9/DisableCopyOnReadDisableCopyOnReadEread_9_disablecopyonread_sparse_model_24_batch_normalization_289_beta"/device:CPU:0*
_output_shapes
 �
Read_9/ReadVariableOpReadVariableOpEread_9_disablecopyonread_sparse_model_24_batch_normalization_289_beta^Read_9/DisableCopyOnRead"/device:CPU:0*
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
Read_10/DisableCopyOnReadDisableCopyOnReadMread_10_disablecopyonread_sparse_model_24_batch_normalization_289_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_10/ReadVariableOpReadVariableOpMread_10_disablecopyonread_sparse_model_24_batch_normalization_289_moving_mean^Read_10/DisableCopyOnRead"/device:CPU:0*
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
Read_11/DisableCopyOnReadDisableCopyOnReadQread_11_disablecopyonread_sparse_model_24_batch_normalization_289_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_11/ReadVariableOpReadVariableOpQread_11_disablecopyonread_sparse_model_24_batch_normalization_289_moving_variance^Read_11/DisableCopyOnRead"/device:CPU:0*
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
Read_12/DisableCopyOnReadDisableCopyOnRead:read_12_disablecopyonread_sparse_model_24_dense_422_kernel"/device:CPU:0*
_output_shapes
 �
Read_12/ReadVariableOpReadVariableOp:read_12_disablecopyonread_sparse_model_24_dense_422_kernel^Read_12/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_24IdentityRead_12/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_25IdentityIdentity_24:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_13/DisableCopyOnReadDisableCopyOnRead8read_13_disablecopyonread_sparse_model_24_dense_422_bias"/device:CPU:0*
_output_shapes
 �
Read_13/ReadVariableOpReadVariableOp8read_13_disablecopyonread_sparse_model_24_dense_422_bias^Read_13/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_26IdentityRead_13/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_27IdentityIdentity_26:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_14/DisableCopyOnReadDisableCopyOnReadGread_14_disablecopyonread_sparse_model_24_batch_normalization_290_gamma"/device:CPU:0*
_output_shapes
 �
Read_14/ReadVariableOpReadVariableOpGread_14_disablecopyonread_sparse_model_24_batch_normalization_290_gamma^Read_14/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_28IdentityRead_14/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_29IdentityIdentity_28:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_15/DisableCopyOnReadDisableCopyOnReadFread_15_disablecopyonread_sparse_model_24_batch_normalization_290_beta"/device:CPU:0*
_output_shapes
 �
Read_15/ReadVariableOpReadVariableOpFread_15_disablecopyonread_sparse_model_24_batch_normalization_290_beta^Read_15/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_30IdentityRead_15/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_31IdentityIdentity_30:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_16/DisableCopyOnReadDisableCopyOnReadMread_16_disablecopyonread_sparse_model_24_batch_normalization_290_moving_mean"/device:CPU:0*
_output_shapes
 �
Read_16/ReadVariableOpReadVariableOpMread_16_disablecopyonread_sparse_model_24_batch_normalization_290_moving_mean^Read_16/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_32IdentityRead_16/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_33IdentityIdentity_32:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_17/DisableCopyOnReadDisableCopyOnReadQread_17_disablecopyonread_sparse_model_24_batch_normalization_290_moving_variance"/device:CPU:0*
_output_shapes
 �
Read_17/ReadVariableOpReadVariableOpQread_17_disablecopyonread_sparse_model_24_batch_normalization_290_moving_variance^Read_17/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_34IdentityRead_17/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_35IdentityIdentity_34:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_18/DisableCopyOnReadDisableCopyOnRead:read_18_disablecopyonread_sparse_model_24_dense_423_kernel"/device:CPU:0*
_output_shapes
 �
Read_18/ReadVariableOpReadVariableOp:read_18_disablecopyonread_sparse_model_24_dense_423_kernel^Read_18/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_36IdentityRead_18/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_37IdentityIdentity_36:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_19/DisableCopyOnReadDisableCopyOnRead8read_19_disablecopyonread_sparse_model_24_dense_423_bias"/device:CPU:0*
_output_shapes
 �
Read_19/ReadVariableOpReadVariableOp8read_19_disablecopyonread_sparse_model_24_dense_423_bias^Read_19/DisableCopyOnRead"/device:CPU:0*
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
Read_22/DisableCopyOnReadDisableCopyOnReadAread_22_disablecopyonread_adam_m_sparse_model_24_dense_420_kernel"/device:CPU:0*
_output_shapes
 �
Read_22/ReadVariableOpReadVariableOpAread_22_disablecopyonread_adam_m_sparse_model_24_dense_420_kernel^Read_22/DisableCopyOnRead"/device:CPU:0*
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
Read_23/DisableCopyOnReadDisableCopyOnReadAread_23_disablecopyonread_adam_v_sparse_model_24_dense_420_kernel"/device:CPU:0*
_output_shapes
 �
Read_23/ReadVariableOpReadVariableOpAread_23_disablecopyonread_adam_v_sparse_model_24_dense_420_kernel^Read_23/DisableCopyOnRead"/device:CPU:0*
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
Read_24/DisableCopyOnReadDisableCopyOnRead?read_24_disablecopyonread_adam_m_sparse_model_24_dense_420_bias"/device:CPU:0*
_output_shapes
 �
Read_24/ReadVariableOpReadVariableOp?read_24_disablecopyonread_adam_m_sparse_model_24_dense_420_bias^Read_24/DisableCopyOnRead"/device:CPU:0*
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
Read_25/DisableCopyOnReadDisableCopyOnRead?read_25_disablecopyonread_adam_v_sparse_model_24_dense_420_bias"/device:CPU:0*
_output_shapes
 �
Read_25/ReadVariableOpReadVariableOp?read_25_disablecopyonread_adam_v_sparse_model_24_dense_420_bias^Read_25/DisableCopyOnRead"/device:CPU:0*
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
Read_26/DisableCopyOnReadDisableCopyOnReadNread_26_disablecopyonread_adam_m_sparse_model_24_batch_normalization_288_gamma"/device:CPU:0*
_output_shapes
 �
Read_26/ReadVariableOpReadVariableOpNread_26_disablecopyonread_adam_m_sparse_model_24_batch_normalization_288_gamma^Read_26/DisableCopyOnRead"/device:CPU:0*
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
Read_27/DisableCopyOnReadDisableCopyOnReadNread_27_disablecopyonread_adam_v_sparse_model_24_batch_normalization_288_gamma"/device:CPU:0*
_output_shapes
 �
Read_27/ReadVariableOpReadVariableOpNread_27_disablecopyonread_adam_v_sparse_model_24_batch_normalization_288_gamma^Read_27/DisableCopyOnRead"/device:CPU:0*
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
Read_28/DisableCopyOnReadDisableCopyOnReadMread_28_disablecopyonread_adam_m_sparse_model_24_batch_normalization_288_beta"/device:CPU:0*
_output_shapes
 �
Read_28/ReadVariableOpReadVariableOpMread_28_disablecopyonread_adam_m_sparse_model_24_batch_normalization_288_beta^Read_28/DisableCopyOnRead"/device:CPU:0*
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
Read_29/DisableCopyOnReadDisableCopyOnReadMread_29_disablecopyonread_adam_v_sparse_model_24_batch_normalization_288_beta"/device:CPU:0*
_output_shapes
 �
Read_29/ReadVariableOpReadVariableOpMread_29_disablecopyonread_adam_v_sparse_model_24_batch_normalization_288_beta^Read_29/DisableCopyOnRead"/device:CPU:0*
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
Read_30/DisableCopyOnReadDisableCopyOnReadAread_30_disablecopyonread_adam_m_sparse_model_24_dense_421_kernel"/device:CPU:0*
_output_shapes
 �
Read_30/ReadVariableOpReadVariableOpAread_30_disablecopyonread_adam_m_sparse_model_24_dense_421_kernel^Read_30/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_60IdentityRead_30/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_61IdentityIdentity_60:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_31/DisableCopyOnReadDisableCopyOnReadAread_31_disablecopyonread_adam_v_sparse_model_24_dense_421_kernel"/device:CPU:0*
_output_shapes
 �
Read_31/ReadVariableOpReadVariableOpAread_31_disablecopyonread_adam_v_sparse_model_24_dense_421_kernel^Read_31/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_62IdentityRead_31/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_63IdentityIdentity_62:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_32/DisableCopyOnReadDisableCopyOnRead?read_32_disablecopyonread_adam_m_sparse_model_24_dense_421_bias"/device:CPU:0*
_output_shapes
 �
Read_32/ReadVariableOpReadVariableOp?read_32_disablecopyonread_adam_m_sparse_model_24_dense_421_bias^Read_32/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_64IdentityRead_32/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_65IdentityIdentity_64:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_33/DisableCopyOnReadDisableCopyOnRead?read_33_disablecopyonread_adam_v_sparse_model_24_dense_421_bias"/device:CPU:0*
_output_shapes
 �
Read_33/ReadVariableOpReadVariableOp?read_33_disablecopyonread_adam_v_sparse_model_24_dense_421_bias^Read_33/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_66IdentityRead_33/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_67IdentityIdentity_66:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_34/DisableCopyOnReadDisableCopyOnReadNread_34_disablecopyonread_adam_m_sparse_model_24_batch_normalization_289_gamma"/device:CPU:0*
_output_shapes
 �
Read_34/ReadVariableOpReadVariableOpNread_34_disablecopyonread_adam_m_sparse_model_24_batch_normalization_289_gamma^Read_34/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_68IdentityRead_34/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_69IdentityIdentity_68:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_35/DisableCopyOnReadDisableCopyOnReadNread_35_disablecopyonread_adam_v_sparse_model_24_batch_normalization_289_gamma"/device:CPU:0*
_output_shapes
 �
Read_35/ReadVariableOpReadVariableOpNread_35_disablecopyonread_adam_v_sparse_model_24_batch_normalization_289_gamma^Read_35/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_70IdentityRead_35/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_71IdentityIdentity_70:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_36/DisableCopyOnReadDisableCopyOnReadMread_36_disablecopyonread_adam_m_sparse_model_24_batch_normalization_289_beta"/device:CPU:0*
_output_shapes
 �
Read_36/ReadVariableOpReadVariableOpMread_36_disablecopyonread_adam_m_sparse_model_24_batch_normalization_289_beta^Read_36/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_72IdentityRead_36/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_73IdentityIdentity_72:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_37/DisableCopyOnReadDisableCopyOnReadMread_37_disablecopyonread_adam_v_sparse_model_24_batch_normalization_289_beta"/device:CPU:0*
_output_shapes
 �
Read_37/ReadVariableOpReadVariableOpMread_37_disablecopyonread_adam_v_sparse_model_24_batch_normalization_289_beta^Read_37/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_74IdentityRead_37/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_75IdentityIdentity_74:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_38/DisableCopyOnReadDisableCopyOnReadAread_38_disablecopyonread_adam_m_sparse_model_24_dense_422_kernel"/device:CPU:0*
_output_shapes
 �
Read_38/ReadVariableOpReadVariableOpAread_38_disablecopyonread_adam_m_sparse_model_24_dense_422_kernel^Read_38/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_76IdentityRead_38/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_77IdentityIdentity_76:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_39/DisableCopyOnReadDisableCopyOnReadAread_39_disablecopyonread_adam_v_sparse_model_24_dense_422_kernel"/device:CPU:0*
_output_shapes
 �
Read_39/ReadVariableOpReadVariableOpAread_39_disablecopyonread_adam_v_sparse_model_24_dense_422_kernel^Read_39/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_78IdentityRead_39/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_79IdentityIdentity_78:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_40/DisableCopyOnReadDisableCopyOnRead?read_40_disablecopyonread_adam_m_sparse_model_24_dense_422_bias"/device:CPU:0*
_output_shapes
 �
Read_40/ReadVariableOpReadVariableOp?read_40_disablecopyonread_adam_m_sparse_model_24_dense_422_bias^Read_40/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_80IdentityRead_40/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_81IdentityIdentity_80:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_41/DisableCopyOnReadDisableCopyOnRead?read_41_disablecopyonread_adam_v_sparse_model_24_dense_422_bias"/device:CPU:0*
_output_shapes
 �
Read_41/ReadVariableOpReadVariableOp?read_41_disablecopyonread_adam_v_sparse_model_24_dense_422_bias^Read_41/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_82IdentityRead_41/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_83IdentityIdentity_82:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_42/DisableCopyOnReadDisableCopyOnReadNread_42_disablecopyonread_adam_m_sparse_model_24_batch_normalization_290_gamma"/device:CPU:0*
_output_shapes
 �
Read_42/ReadVariableOpReadVariableOpNread_42_disablecopyonread_adam_m_sparse_model_24_batch_normalization_290_gamma^Read_42/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_84IdentityRead_42/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_85IdentityIdentity_84:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_43/DisableCopyOnReadDisableCopyOnReadNread_43_disablecopyonread_adam_v_sparse_model_24_batch_normalization_290_gamma"/device:CPU:0*
_output_shapes
 �
Read_43/ReadVariableOpReadVariableOpNread_43_disablecopyonread_adam_v_sparse_model_24_batch_normalization_290_gamma^Read_43/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_86IdentityRead_43/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_87IdentityIdentity_86:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_44/DisableCopyOnReadDisableCopyOnReadMread_44_disablecopyonread_adam_m_sparse_model_24_batch_normalization_290_beta"/device:CPU:0*
_output_shapes
 �
Read_44/ReadVariableOpReadVariableOpMread_44_disablecopyonread_adam_m_sparse_model_24_batch_normalization_290_beta^Read_44/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_88IdentityRead_44/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_89IdentityIdentity_88:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_45/DisableCopyOnReadDisableCopyOnReadMread_45_disablecopyonread_adam_v_sparse_model_24_batch_normalization_290_beta"/device:CPU:0*
_output_shapes
 �
Read_45/ReadVariableOpReadVariableOpMread_45_disablecopyonread_adam_v_sparse_model_24_batch_normalization_290_beta^Read_45/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:*
dtype0k
Identity_90IdentityRead_45/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:a
Identity_91IdentityIdentity_90:output:0"/device:CPU:0*
T0*
_output_shapes
:�
Read_46/DisableCopyOnReadDisableCopyOnReadAread_46_disablecopyonread_adam_m_sparse_model_24_dense_423_kernel"/device:CPU:0*
_output_shapes
 �
Read_46/ReadVariableOpReadVariableOpAread_46_disablecopyonread_adam_m_sparse_model_24_dense_423_kernel^Read_46/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_92IdentityRead_46/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_93IdentityIdentity_92:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_47/DisableCopyOnReadDisableCopyOnReadAread_47_disablecopyonread_adam_v_sparse_model_24_dense_423_kernel"/device:CPU:0*
_output_shapes
 �
Read_47/ReadVariableOpReadVariableOpAread_47_disablecopyonread_adam_v_sparse_model_24_dense_423_kernel^Read_47/DisableCopyOnRead"/device:CPU:0*
_output_shapes

:*
dtype0o
Identity_94IdentityRead_47/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

:e
Identity_95IdentityIdentity_94:output:0"/device:CPU:0*
T0*
_output_shapes

:�
Read_48/DisableCopyOnReadDisableCopyOnRead?read_48_disablecopyonread_adam_m_sparse_model_24_dense_423_bias"/device:CPU:0*
_output_shapes
 �
Read_48/ReadVariableOpReadVariableOp?read_48_disablecopyonread_adam_m_sparse_model_24_dense_423_bias^Read_48/DisableCopyOnRead"/device:CPU:0*
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
Read_49/DisableCopyOnReadDisableCopyOnRead?read_49_disablecopyonread_adam_v_sparse_model_24_dense_423_bias"/device:CPU:0*
_output_shapes
 �
Read_49/ReadVariableOpReadVariableOp?read_49_disablecopyonread_adam_v_sparse_model_24_dense_423_bias^Read_49/DisableCopyOnRead"/device:CPU:0*
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
_user_specified_name	total_1:E2A
?
_user_specified_name'%Adam/v/sparse_model_24/dense_423/bias:E1A
?
_user_specified_name'%Adam/m/sparse_model_24/dense_423/bias:G0C
A
_user_specified_name)'Adam/v/sparse_model_24/dense_423/kernel:G/C
A
_user_specified_name)'Adam/m/sparse_model_24/dense_423/kernel:S.O
M
_user_specified_name53Adam/v/sparse_model_24/batch_normalization_290/beta:S-O
M
_user_specified_name53Adam/m/sparse_model_24/batch_normalization_290/beta:T,P
N
_user_specified_name64Adam/v/sparse_model_24/batch_normalization_290/gamma:T+P
N
_user_specified_name64Adam/m/sparse_model_24/batch_normalization_290/gamma:E*A
?
_user_specified_name'%Adam/v/sparse_model_24/dense_422/bias:E)A
?
_user_specified_name'%Adam/m/sparse_model_24/dense_422/bias:G(C
A
_user_specified_name)'Adam/v/sparse_model_24/dense_422/kernel:G'C
A
_user_specified_name)'Adam/m/sparse_model_24/dense_422/kernel:S&O
M
_user_specified_name53Adam/v/sparse_model_24/batch_normalization_289/beta:S%O
M
_user_specified_name53Adam/m/sparse_model_24/batch_normalization_289/beta:T$P
N
_user_specified_name64Adam/v/sparse_model_24/batch_normalization_289/gamma:T#P
N
_user_specified_name64Adam/m/sparse_model_24/batch_normalization_289/gamma:E"A
?
_user_specified_name'%Adam/v/sparse_model_24/dense_421/bias:E!A
?
_user_specified_name'%Adam/m/sparse_model_24/dense_421/bias:G C
A
_user_specified_name)'Adam/v/sparse_model_24/dense_421/kernel:GC
A
_user_specified_name)'Adam/m/sparse_model_24/dense_421/kernel:SO
M
_user_specified_name53Adam/v/sparse_model_24/batch_normalization_288/beta:SO
M
_user_specified_name53Adam/m/sparse_model_24/batch_normalization_288/beta:TP
N
_user_specified_name64Adam/v/sparse_model_24/batch_normalization_288/gamma:TP
N
_user_specified_name64Adam/m/sparse_model_24/batch_normalization_288/gamma:EA
?
_user_specified_name'%Adam/v/sparse_model_24/dense_420/bias:EA
?
_user_specified_name'%Adam/m/sparse_model_24/dense_420/bias:GC
A
_user_specified_name)'Adam/v/sparse_model_24/dense_420/kernel:GC
A
_user_specified_name)'Adam/m/sparse_model_24/dense_420/kernel:-)
'
_user_specified_namelearning_rate:)%
#
_user_specified_name	iteration:>:
8
_user_specified_name sparse_model_24/dense_423/bias:@<
:
_user_specified_name" sparse_model_24/dense_423/kernel:WS
Q
_user_specified_name97sparse_model_24/batch_normalization_290/moving_variance:SO
M
_user_specified_name53sparse_model_24/batch_normalization_290/moving_mean:LH
F
_user_specified_name.,sparse_model_24/batch_normalization_290/beta:MI
G
_user_specified_name/-sparse_model_24/batch_normalization_290/gamma:>:
8
_user_specified_name sparse_model_24/dense_422/bias:@<
:
_user_specified_name" sparse_model_24/dense_422/kernel:WS
Q
_user_specified_name97sparse_model_24/batch_normalization_289/moving_variance:SO
M
_user_specified_name53sparse_model_24/batch_normalization_289/moving_mean:L
H
F
_user_specified_name.,sparse_model_24/batch_normalization_289/beta:M	I
G
_user_specified_name/-sparse_model_24/batch_normalization_289/gamma:>:
8
_user_specified_name sparse_model_24/dense_421/bias:@<
:
_user_specified_name" sparse_model_24/dense_421/kernel:WS
Q
_user_specified_name97sparse_model_24/batch_normalization_288/moving_variance:SO
M
_user_specified_name53sparse_model_24/batch_normalization_288/moving_mean:LH
F
_user_specified_name.,sparse_model_24/batch_normalization_288/beta:MI
G
_user_specified_name/-sparse_model_24/batch_normalization_288/gamma:>:
8
_user_specified_name sparse_model_24/dense_420/bias:@<
:
_user_specified_name" sparse_model_24/dense_420/kernel:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix"�L
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
0__inference_sparse_model_24_layer_call_fn_525808
0__inference_sparse_model_24_layer_call_fn_525853�
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
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525712
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525763�
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
!__inference__wrapped_model_525380input_1"�
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
3:1	�2 sparse_model_24/dense_420/kernel
,:*2sparse_model_24/dense_420/bias
;:92-sparse_model_24/batch_normalization_288/gamma
::82,sparse_model_24/batch_normalization_288/beta
C:A (23sparse_model_24/batch_normalization_288/moving_mean
G:E (27sparse_model_24/batch_normalization_288/moving_variance
2:02 sparse_model_24/dense_421/kernel
,:*2sparse_model_24/dense_421/bias
;:92-sparse_model_24/batch_normalization_289/gamma
::82,sparse_model_24/batch_normalization_289/beta
C:A (23sparse_model_24/batch_normalization_289/moving_mean
G:E (27sparse_model_24/batch_normalization_289/moving_variance
2:02 sparse_model_24/dense_422/kernel
,:*2sparse_model_24/dense_422/bias
;:92-sparse_model_24/batch_normalization_290/gamma
::82,sparse_model_24/batch_normalization_290/beta
C:A (23sparse_model_24/batch_normalization_290/moving_mean
G:E (27sparse_model_24/batch_normalization_290/moving_variance
2:02 sparse_model_24/dense_423/kernel
,:*2sparse_model_24/dense_423/bias
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
0__inference_sparse_model_24_layer_call_fn_525808input_1"�
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
0__inference_sparse_model_24_layer_call_fn_525853input_1"�
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
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525712input_1"�
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
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525763input_1"�
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
*__inference_dense_420_layer_call_fn_525946�
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
E__inference_dense_420_layer_call_and_return_conditional_losses_525956�
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
8__inference_batch_normalization_288_layer_call_fn_525969
8__inference_batch_normalization_288_layer_call_fn_525982�
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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526016
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526036�
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
*__inference_dense_421_layer_call_fn_526045�
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
E__inference_dense_421_layer_call_and_return_conditional_losses_526055�
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
8__inference_batch_normalization_289_layer_call_fn_526068
8__inference_batch_normalization_289_layer_call_fn_526081�
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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526115
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526135�
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
*__inference_dense_422_layer_call_fn_526144�
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
E__inference_dense_422_layer_call_and_return_conditional_losses_526154�
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
8__inference_batch_normalization_290_layer_call_fn_526167
8__inference_batch_normalization_290_layer_call_fn_526180�
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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526214
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526234�
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
*__inference_dense_423_layer_call_fn_526243�
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
E__inference_dense_423_layer_call_and_return_conditional_losses_526254�
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
$__inference_signature_wrapper_525937input_1"�
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
*__inference_dense_420_layer_call_fn_525946inputs"�
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
E__inference_dense_420_layer_call_and_return_conditional_losses_525956inputs"�
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
8__inference_batch_normalization_288_layer_call_fn_525969inputs"�
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
8__inference_batch_normalization_288_layer_call_fn_525982inputs"�
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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526016inputs"�
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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526036inputs"�
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
*__inference_dense_421_layer_call_fn_526045inputs"�
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
E__inference_dense_421_layer_call_and_return_conditional_losses_526055inputs"�
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
8__inference_batch_normalization_289_layer_call_fn_526068inputs"�
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
8__inference_batch_normalization_289_layer_call_fn_526081inputs"�
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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526115inputs"�
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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526135inputs"�
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
*__inference_dense_422_layer_call_fn_526144inputs"�
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
E__inference_dense_422_layer_call_and_return_conditional_losses_526154inputs"�
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
8__inference_batch_normalization_290_layer_call_fn_526167inputs"�
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
8__inference_batch_normalization_290_layer_call_fn_526180inputs"�
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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526214inputs"�
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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526234inputs"�
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
*__inference_dense_423_layer_call_fn_526243inputs"�
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
E__inference_dense_423_layer_call_and_return_conditional_losses_526254inputs"�
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
8:6	�2'Adam/m/sparse_model_24/dense_420/kernel
8:6	�2'Adam/v/sparse_model_24/dense_420/kernel
1:/2%Adam/m/sparse_model_24/dense_420/bias
1:/2%Adam/v/sparse_model_24/dense_420/bias
@:>24Adam/m/sparse_model_24/batch_normalization_288/gamma
@:>24Adam/v/sparse_model_24/batch_normalization_288/gamma
?:=23Adam/m/sparse_model_24/batch_normalization_288/beta
?:=23Adam/v/sparse_model_24/batch_normalization_288/beta
7:52'Adam/m/sparse_model_24/dense_421/kernel
7:52'Adam/v/sparse_model_24/dense_421/kernel
1:/2%Adam/m/sparse_model_24/dense_421/bias
1:/2%Adam/v/sparse_model_24/dense_421/bias
@:>24Adam/m/sparse_model_24/batch_normalization_289/gamma
@:>24Adam/v/sparse_model_24/batch_normalization_289/gamma
?:=23Adam/m/sparse_model_24/batch_normalization_289/beta
?:=23Adam/v/sparse_model_24/batch_normalization_289/beta
7:52'Adam/m/sparse_model_24/dense_422/kernel
7:52'Adam/v/sparse_model_24/dense_422/kernel
1:/2%Adam/m/sparse_model_24/dense_422/bias
1:/2%Adam/v/sparse_model_24/dense_422/bias
@:>24Adam/m/sparse_model_24/batch_normalization_290/gamma
@:>24Adam/v/sparse_model_24/batch_normalization_290/gamma
?:=23Adam/m/sparse_model_24/batch_normalization_290/beta
?:=23Adam/v/sparse_model_24/batch_normalization_290/beta
7:52'Adam/m/sparse_model_24/dense_423/kernel
7:52'Adam/v/sparse_model_24/dense_423/kernel
1:/2%Adam/m/sparse_model_24/dense_423/bias
1:/2%Adam/v/sparse_model_24/dense_423/bias
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
!__inference__wrapped_model_525380~!" #$1�.
'�$
"�
input_1����������
� "3�0
.
output_1"�
output_1����������
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526016m7�4
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
S__inference_batch_normalization_288_layer_call_and_return_conditional_losses_526036m7�4
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
8__inference_batch_normalization_288_layer_call_fn_525969b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_288_layer_call_fn_525982b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526115m7�4
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
S__inference_batch_normalization_289_layer_call_and_return_conditional_losses_526135m7�4
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
8__inference_batch_normalization_289_layer_call_fn_526068b7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_289_layer_call_fn_526081b7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526214m!" 7�4
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
S__inference_batch_normalization_290_layer_call_and_return_conditional_losses_526234m!" 7�4
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
8__inference_batch_normalization_290_layer_call_fn_526167b!" 7�4
-�*
 �
inputs���������
p

 
� "!�
unknown����������
8__inference_batch_normalization_290_layer_call_fn_526180b!" 7�4
-�*
 �
inputs���������
p 

 
� "!�
unknown����������
E__inference_dense_420_layer_call_and_return_conditional_losses_525956d0�-
&�#
!�
inputs����������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_420_layer_call_fn_525946Y0�-
&�#
!�
inputs����������
� "!�
unknown����������
E__inference_dense_421_layer_call_and_return_conditional_losses_526055c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_421_layer_call_fn_526045X/�,
%�"
 �
inputs���������
� "!�
unknown����������
E__inference_dense_422_layer_call_and_return_conditional_losses_526154c/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_422_layer_call_fn_526144X/�,
%�"
 �
inputs���������
� "!�
unknown����������
E__inference_dense_423_layer_call_and_return_conditional_losses_526254c#$/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������
� �
*__inference_dense_423_layer_call_fn_526243X#$/�,
%�"
 �
inputs���������
� "!�
unknown����������
$__inference_signature_wrapper_525937�!" #$<�9
� 
2�/
-
input_1"�
input_1����������"3�0
.
output_1"�
output_1����������
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525712�!" #$A�>
'�$
"�
input_1����������
�

trainingp",�)
"�
tensor_0���������
� �
K__inference_sparse_model_24_layer_call_and_return_conditional_losses_525763�!" #$A�>
'�$
"�
input_1����������
�

trainingp ",�)
"�
tensor_0���������
� �
0__inference_sparse_model_24_layer_call_fn_525808|!" #$A�>
'�$
"�
input_1����������
�

trainingp"!�
unknown����������
0__inference_sparse_model_24_layer_call_fn_525853|!" #$A�>
'�$
"�
input_1����������
�

trainingp "!�
unknown���������