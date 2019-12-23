# Principal Component Analysis using Linear Algebra - Step by Step Guide

I was using PCA from quite a long time but a few days before I tried to learn how it works. So in this blog I have tried to explain PCA in an intuitive way, you will be able to understand the mathematics and will also see how and why do we use it.

Before dive into PCA we need to understand SVD, Eigenvalues and Eigenvectors.

## what is **Single Value Decomposition.**

Lets start with one example.
$$\vec X = \begin{bmatrix}1 \\ 3\end{bmatrix} A = \begin{bmatrix}2 & 1\\-1 & 1\end{bmatrix}$$
$$Y = A\vec X = \begin{bmatrix}2 & 1\\-1 & 1\end{bmatrix} \begin{bmatrix}1 \\ 3\end{bmatrix} = \begin{bmatrix}5 \\ 2\end{bmatrix}$$

We haven't done anything special but let we try to see graphically.
* Graphically we can say we transformed one vector to another (rotated and stratched). This is what $A$ does.
* $A = \begin{bmatrix}\cos\theta & -\sin\theta\\\sin\theta & \cos\theta\end{bmatrix} \leftarrow$ rotation matrix, if you want to rotate any vector to angle $\theta$
*  $A = \begin{bmatrix}\alpha & 0\\0 & \alpha  \end{bmatrix} \leftarrow$ $\alpha$ is stretching.

Now we can think that $A$ is made up of some combination of rotation & stretching.

Let's see only stretching:
$$ = \begin{bmatrix}2 & 0\\0 & 2  \end{bmatrix} \begin{bmatrix}1\\3  \end{bmatrix} =  \begin{bmatrix}2\\6  \end{bmatrix}$$

When we say of a single vector it is represented by a line. Now Let say there are 2-Dimensions(orthonormal) representing a circle.

$A$ is matrix represented as : $A \in \mathbb{C}^{m\times n}$

Multiplication with metrix $A$

* of vector(1-D) rotates and stretches it but what will happen to Circle (2-D), Sphere (3-D), hyper Sphere(>3-D).

* So, Circle will be converted into Ellipse, Sphere into Ellipsoid and Hyper-sphere into hyper-ellipsoid/ellipse.

For circle we will get an ellipse, with Major and Minor axes, we will be having a new co-ordinate system, it will be represented by ($\vec u\sigma$)

Unit vectors ($\vec u_1, \vec u_2$) tells us, in which direction you are going & $\sigma_1, \sigma_2$ tells, how much we are going in that perticular direction.

We started with one vector space and after multipling it with $A$ we ended up in new vector space.
$$\vec v_1, \vec v_2, ... , \vec v_n \rightarrow \vec u_1\sigma_1, \vec u_2\sigma_2, ... , \vec u_n\sigma_n$$

For this operation our matrix $A$ should be rank efficient, means members of rows and columns are linearly independent.

 $\vec u_1\sigma_1, \vec u_2\sigma_2$ are principal axes and stretching factors $\sigma_1, \sigma_2$ are singular values.
$$A \vec v_1 = \sigma_1 \vec u_1$$
We took a matrix $A$(rotation, stretch) on vector $\vec v_1$ and produced $\sigma_1 \vec u_1$. This is similar to eigenvalue problem with a little difference.

Eigenvalue problem: $A \vec v = \lambda \vec v$ (both sides $\vec v$)

Try to generalise:
$$A \vec v_j = \sigma_j\vec u_j \text{ where j = 1, 2, 3, 4,..., n}$$

$$\begin{bmatrix} \\ \\A\\ \\ \\ \end{bmatrix}\begin{bmatrix}\vec v_1 & \vec v_2 & ... & \vec v_n\end{bmatrix} = \begin{bmatrix} \\ \vec u_1 & \vec u_2 & ... & \vec u_n \\ \\\end{bmatrix}\begin{bmatrix} \sigma_1 & 0 & ... & 0 \\ 0 & \sigma_2 & ... & 0 \\ ... & ... & ... & ... \\ 0 & 0 & ... & \sigma_n \end{bmatrix}$$
$$A_{m\times n}V_{n\times n} = \hat U_{Rotation} \hat \Sigma_{Stretching}$$

From one orthonormal basis to another orthonormal basis.

All stretching is picked up by $\hat \Sigma$ and all rotation picked up by $\hat U$, which was produced by $A$.

* Rotation matrix is also called **unitary transformation**.

* And if these are orthogonal components then:
    * inverse of the matrix is same as its transpose.
    * $U^{-1} = U^{*}$
    * same for $V^{-1} = V^{*}$

$$AV = \hat U \hat \Sigma$$
Multiply both sides by $V^{-1}$
$$AVV^{-1} = \hat U \hat \Sigma V^{-1}$$
$$A = \hat U \hat \Sigma V^{*}$$
This is **reduced singular value decomposition.**
$$A_{m \times n} = \hat U_{m \times n} \hat \Sigma_{n \times n} V_{n \times n}^*$$
Now add ($m-n$) silent column to $\hat U$ to make it $(m \times m)$ and ($m-n$) silent rows to $\hat \Sigma$ to make it $(m \times n)$

Also remember diagonal elements of $\Sigma$ are as follows:

$\sigma_1 \geq \sigma_2 \geq ... \geq\sigma_n$

Theorem(SVD): Every matrix $A \in \mathbb{C}^{m \times n}$ has a singular value decomposition.
* Singular values {$\sigma_j$} are uniquely determined, and if $A$ is square then $\sigma_j$ are distinct.
* Singular vectors {$\vec u_j$} and {$\vec v_j$} are also unique up to a complex sign.

Now, how to compute all the things in $A = U \Sigma V^*$ ?

Let's start with $A^TA$
$$A^TA = (U \Sigma V^*)^T(U \Sigma V^*)$$
$$= V \Sigma U^* U \Sigma V^*$$
$$= V \Sigma I \Sigma V^*$$
$$= V \Sigma^2 V^*$$
Now $A^TAV= V \Sigma^2 V^*V$
$$A^TAV= V \Sigma^2 $$
This is same as Eigen value problem $A \vec X = \lambda \vec X$

If we solve eigen value problem we will get $\lambda_j = \sigma_j^2$

How to get $V$ ?
$$AA^T = (U \Sigma V^*)(U \Sigma V^*)^T$$
$$ = U \Sigma V^* V \Sigma U^*$$
$$ = U \Sigma ^2 U^*$$
$AA^TU = U \Sigma ^2 U^*U$
$$AA^TU = U \Sigma ^2$$
Eigen value problem. $U$ is eigen values. Its called self-adjoint or Hermition.

## Eigenvalues and Eigenvectors

The word **Eigen** means **latent** or **characteristic**.
Lets say 12 can be represented by $2\times2\times3$. Similarly a matrix is also can be represented by same vectors, these vectors are called eigen vectors.

Let's see more intuitively:

We have an equation:
$$Ax = \lambda x$$
It is true for **special** x(vectors) and **special** $\lambda$ (numbers).

Let's begin with some numerical example:
$$A = \begin{bmatrix}3 & -1 \\ -1 & 3\end{bmatrix}$$
and try $x =\begin{bmatrix}1 \\ 0\end{bmatrix}$:
$$Ax = \begin{bmatrix}3 \\ -1\end{bmatrix}$$
So, for the above $x$, $x$ and $Ax$ are not in same direction.
Lets try with $x =\begin{bmatrix}0 \\ 1\end{bmatrix}$
$$Ax = \begin{bmatrix}-1 \\ 3\end{bmatrix}$$
Also, for the above $x$, $x$ and $Ax$ are not in same direction.

We can say no $\lambda$ is going to satisfy the given equation for given values.

Now let us take **special** $x =\begin{bmatrix}1 \\ 1\end{bmatrix}$
$$Ax = \begin{bmatrix}2 \\ 2\end{bmatrix}$$
Now both are pointing in same direction, and now we have some value of $\lambda$ that satisfies the given equation: $2\vec x = A\vec x$

We can say this value of $\vec x$ is eigen vector of $A$ and $\lambda =2$ is eigen value.

Take one more **special** $x =\begin{bmatrix}1 \\ -1\end{bmatrix}$
$$Ax = \begin{bmatrix}4 \\ -4\end{bmatrix}$$
Both $Ax$ and $x$ are pointing in same direction. In this case $\lambda =4$. $4\vec x = A\vec x$

So, how to get these special values?
We know $Ax = \lambda x$ is same as $Ax = \lambda I x$

$(A - \lambda I)x = 0$

For above equation

case1:  $x=0$ (not-interesting)

case2: $x\neq0$, now $\det(A-\lambda I) =0$
 * When $\det(A-\lambda I) =0$, $(A-\lambda I)$  is **singular** matrix.
 * $\det(A-\lambda I) =0$ is called **characteristic** equation, and its a polynomial equation and roots $\lambda$ are eigenvalues.

$\det(A-\lambda I) =0$
$$A - \lambda = \begin{bmatrix}3 & -1 \\ -1 & 3\end{bmatrix} - \begin{bmatrix}\lambda & 0 \\ 0 & \lambda\end{bmatrix}$$
$$ = \begin{bmatrix}3 -\lambda & -1 \\ -1 & 3-\lambda\end{bmatrix}$$

Step 1: Compute $\lambda$
$$\det(A-\lambda I) = (3-\lambda)(3-\lambda) - 1$$
$$ = \lambda^2-6\lambda + 8 = 0$$
$$ = (\lambda-4)(\lambda-2) = 0$$
$\lambda =4,2$ are solutions to the equations and are eigen values.

Step 2: Compute eigen vector given eigen values.
Solve for $\vec x$:
For $\lambda=2$ :
$$A -2I = \begin{bmatrix}3-2 & -1 \\-1 & 3-2\end{bmatrix} =  \begin{bmatrix}1 & -1 \\-1 & 1\end{bmatrix}$$
$$\begin{bmatrix}1 & -1 \\-1 & 1\end{bmatrix}\begin{bmatrix} x_1 \\ x_2\end{bmatrix} = 0$$
$$\begin{cases}x_1 - x_2 = 0\\-x_1 + x_2 = 0\end{cases}$$
$$x_1 = x_2$$
so,  all vectors in direction $\begin{bmatrix}1 \\1\end{bmatrix}$ are eigenvector so we call it eigen direction.

For $\lambda=4$ :
$$A -4I = \begin{bmatrix}3-4 & -1 \\-1 & 3-4\end{bmatrix} =  \begin{bmatrix}1 & -1 \\-1 & 1\end{bmatrix}$$
$$\begin{bmatrix}-1 & -1 \\-1 & -1\end{bmatrix}\begin{bmatrix} x_1 \\ x_2\end{bmatrix} = 0$$
$$\begin{cases}x_1 - x_2 = 0\\x_1 - x_2 = 0\end{cases}$$
$$x_1 = -x_2$$
so,  all vectors in direction $\begin{bmatrix}1 \\-1\end{bmatrix}$ are eigenvector so we call it eigen direction.

$[T, D] = eig(A)$, $T$= Vectors, $D$=Values

These satisfies $A = TDT^{-1}$

$T = \begin{bmatrix}  \\x_{\lambda_{1}} & x_{\lambda_{2}} & ... & x_{\lambda_{n}}\\ \\\end{bmatrix}$
$D = \begin{bmatrix}  \lambda_1 &0 & ... & 0\\0 &\lambda_2 & ... & 0\\... &... & ... & ... \\0 &0 & ... & \lambda_n\end{bmatrix}$

Its always easy for us to compute $D^N$ easily but $A^N$ is computationally expensive.
$$D^N = \begin{bmatrix}  \lambda_1^N &0 & ... & 0\\0 &\lambda_2^N & ... & 0\\... &... & ... & ... \\0 &0 & ... & \lambda_n^N\end{bmatrix}$$

$A^N = (TDT^{-1})(TDT^{-1})...(TDT^{-1})$, n-times

$A^N = (TD(T^{-1}T)D(T^{-1}...T)DT^{-1})$

$A^N = (TD^NT^{-1})$


Now we are ready to understand PCA.
## Principal Component Analysis
So let’s begin with a fancy Physics example.

![springMass](./Images/springMass.jpg)

As you can see in the figure there is mass $m$ hanging to a spring with spring constant $k$. Now assume that we pulled the mass in vertical direction and released. Due to this action it will start motion in 1-Direction(up-down). So now we are interested in that motion in one direction, and as you all know it can be calculated using the following equations.

$$F = ma \tag{1}$$
$$Force = mass * acceleration$$
Also,
$$-kf(t) = ma$$
$$-\frac{k}{m}f(t) = a$$
$$-\omega^{2}f(t) = \frac{\text{d}f(t)}{\text{d}t^{2}}$$
Where,
$\omega =$ Angular Frequency
$f(t) =$ Function of displacement
$$f(t) = Acos(\omega t + \omega_{0})$$
So in this way we can solve the **Spring Mass System**.  Now let say we don't have $F=ma$, we only had some data. Can we get $f=ma$  just from the data?

Now take the above system and install some cameras around it to capture some data. For our example we can take 3 cameras, each one of them recording some set of data. Now our objective is to figure just from the data can we say someting about the system.

Let say we got data from cameras in the below format:
* Camera1 : $(\vec X_{a}, \vec Y_{a})$
* Camera2 : $(\vec X_{b}, \vec Y_{b})$
* Camera3 : $(\vec X_{c}, \vec Y_{c})$

Our data matrix will look something like:
$$X = \begin{bmatrix}\vec X_{a} \\\vec Y_{a}\\\vec X_{b}\\\vec Y_{b}\\\vec X_{c}\\\vec Y_{c} \end{bmatrix}$$

Whenever we get data from these kind of sources there are two kind of problems associated with it:
1. Noise
2. Redundancy

We will focus on redundancy,
- is my $\vec X$ is independent of $\vec Y$
- All three cameras are giving me the same information just from different angles.

In our case we have 6 sets of data but we only needed 1. So our **PCA** is going to tell us where should we place our single camera and record only a single value.

Before dive into PCA we need to understand some related concepts:
* variance
* Covariance
Say we have two collections of data:

    $\vec a=[a_{1}, a_{2}, ... ,a_{n} ]$ and $\vec b=[b_{1}, b_{2}, ... ,b_{n} ]$
**Covariance** says how these two data sets are related. On the other hand **variance** says how large changes are there in the data set itself.

    variance is given by:
    $$\sigma_{a}^{2} = \frac{1}{n-1}\vec a \vec a^{T}  \text{ subject to mean= 0}\tag{2}$$
    If the product  $\vec a \vec a^{T}$ is big, it means lot of variance. If its small means low variance.
    You can assume a $sine$ wave on $x-axis$ is with mean zero and its amplitude($A$) says its variance, bigger value of $A$ says high variance and smaller value says low variance.

    Covariance is given by:
    $$\sigma_{a,b}^{2} = \frac{1}{n-1}\vec a \vec b^{T}  \text{ subject to mean= 0} \tag{3}$$
    $\vec a \vec b^{T}$ tells you how much $\vec a$ & $\vec b$ are in same direction. Maximum value is **1**. Value =  **1**, means they are in same direction($\theta=0^0$),  **0** means they are orthogonal to each other($\theta=90^0$) and value =  **-1** says they are in opposite direction($\theta=180^0$).


Now let's come back to Spring Mass System. Whenever lots of data is available we try to remove redundancy. For example in our dataset we only need *1-degree of freedom* but we have *6-degrees of freedom* of information. We need to convert it back to *1-degree*.

We can generalize our equation from vectors to metrices as:

$$C_x = \frac{1}{1-n}XX^T$$
In our case $C_x$ is $6\times 6$ matrix.
$$C_x =  \begin{bmatrix}\mathbf{\sigma_{X_{a}X_{a}}^2} & \sigma_{X_{a}Y_{a}}^2 & \sigma_{X_{a}X_{b}}^2 & \sigma_{X_{a}Y_{b}}^2 & \sigma_{X_{a}X_{c}}^2 & \sigma_{X_{a}Y_{c}}^2\\\sigma_{Y_{a}X_{a}}^2 & \mathbf{\sigma_{Y_{a}Y_{a}}^2} & \sigma_{Y_{a}X_{b}}^2 & \sigma_{Y_{a}Y_{b}}^2 & \sigma_{Y_{a}X_{c}}^2 & \sigma_{Y_{a}Y_{c}}^2\\\sigma_{X_{b}X_{a}}^2 & \sigma_{X_{b}Y_{a}}^2 & \mathbf{\sigma_{X_{b}X_{b}}^2} & \sigma_{X_{b}Y_{b}}^2 & \sigma_{X_{b}X_{c}}^2 & \sigma_{X_{b}Y_{c}}^2\\\sigma_{Y_{b}X_{a}}^2 & \sigma_{Y_{b}Y_{a}}^2 & \sigma_{Y_{b}X_{b}}^2 & \mathbf{\sigma_{Y_{b}Y_{b}}^2} & \sigma_{Y_{b}X_{c}}^2 & \sigma_{Y_{b}Y_{c}}^2\\\sigma_{X_{c}X_{a}}^2 & \sigma_{X_{c}Y_{a}}^2 & \sigma_{X_{c}X_{b}}^2 & \sigma_{X_{c}Y_{b}}^2 & \mathbf{\sigma_{X_{c}X_{c}}^2} & \sigma_{X_{c}Y_{c}}^2\\\sigma_{Y_{c}X_{a}}^2 & \sigma_{Y_{c}Y_{a}}^2 & \sigma_{Y_{c}X_{b}}^2 & \sigma_{Y_{c}Y_{b}}^2 & \sigma_{Y_{c}X_{c}}^2 & \mathbf{\sigma_{Y_{c}Y_{c}}^2}\end{bmatrix} \tag{4}$$

What we can see in above matrix is:
* Diagonals: $\rightarrow$ variance measures.
* Off-diagonals: $\rightarrow$ Co-variance between all pairs.
    * All off-diagonal elements are Symmetric, these type of metrices are called **symmetric** metrices or **self-adjoint** metrices or **Hermition** metrices.
       * If values are very small: It means they are Statistcally independent.
       * If values are very big: It means they are Statistically dependent $\rightarrow$ means higher redundancy.

If diagonal values are higher (variance is higher) that means there is a lot of stuff is happening and if its smaller then no much is going on.

Diagonal terms with big values are important $\rightarrow$ all system stuff is happening.

Our main aim is to make $X$ diagonal $\rightarrow$ no redundancy $\rightarrow$  covariance should be minimum ($0$).
$$X \Rightarrow\text{MAKE DIAGONAL}$$

Or we can say is we want to diagonalise the system $\rightarrow$ change the basis we are working in so it becomes diagonal.

If all elements are other than diagonal are zero means there is no dependency, they are independent, there is no redundancy.

After diagonalising the system we will be left with diagonal terms, now we can arrange the diagonal terms in descending order.
$$C_x = \begin{bmatrix}biggest & 0 & 0 &0\\0 & ... & 0 &0 \\0 & 0 & ... &0\\0 &0 & 0 &smallest\end{bmatrix}$$

Biggest diagonals are going to tell us the dynamics are strongest.

Now our system looks like Singular Value Decomposition (SVD).

$X$ is data matrix.

$XX^T = S \Lambda S^{-1}$ eigenvalue decomposition.

Where $S$ is matrix of eigenvectors of $XX^T$.

And for orthonormal metrices $S^{-1} = S^T$.
$\Lambda$ is diagonal matrix with eigenvalues of $XX^T$.

We wanted to create a new frame of reference, to make a new set of measurement.

Let's say:

$Y = S^TX$, so $Y$ is our new set.
$$C_Y = \frac{1}{1-n}YY^T$$
$$C_Y = \frac{1}{1-n}(S^TX)(S^TX)^T$$
$$= \frac{1}{1-n}(S^TX)(X^TS)$$
$$= \frac{1}{1-n}(S^TS \Lambda S^{-1}S)$$
$$C_Y = \frac{\Lambda}{1-n} \text{ where $\Lambda$ is a diagonal matrix.}$$

Now for our dataset:

$$\Lambda =  \begin{bmatrix}\sigma & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 \end{bmatrix}$$


So it says, only one direction matters and also says it is a single degree of freedom problem.

So, this solution was using eigenvalues and eigen vectors.

Now we will solve the above equation by SVD.

$X$ is our data matrix  =  $U\Sigma V^*$
$$Y = U^* X$$
$$C_Y = \frac{1}{1-n} YY^T$$
$$C_Y = \frac{1}{1-n} (U^* X)(U^* X)^T$$
$$C_Y = \frac{1}{1-n} (U^* X)(X^TU )$$
$$C_Y = \frac{1}{1-n} (U^* U\Sigma V^*)((U\Sigma V^*)^TU )$$
$$C_Y = \frac{1}{1-n} (\Sigma V^*)(V\Sigma U^*U )$$
$$C_Y = \frac{1}{1-n} \Sigma^2$$
$\Sigma^2$ is same as $\Lambda \rightarrow$ $\lambda_i = \sigma_i^2$ (Connection between Single value and eigen vectors.)
