Hey Costas,

Glad to see someone using the kmeans stuff.

> --I confess to not understanding the docstring

Sorry for the confusion.  I'll try to explain thing more clearly.  If it works, will use this as the doc. :)

> However, I am not quite following the kmeans() functionality (I am new to 
> this clustering business, so this maybe a stupid newbie question): my docs 
> tell me that kmeans should partition a dataset into k clusters.  So, I 
> expect vq.kmeans(dataset, 2) to return to me dataset split up into two 
> "equivalent" datasets.

Splitting the data into two data sets is actually a two or three step process. 
Here's a complete example.

"Observations" are just another name for a data point.  The obs matrix is a 2D 
array of data points.  For example if our data set includes height, weight, and 
40-yard speed of football players, you might have the following (fictitous) 
data:

obs:
                            lb        inches     seconds       
                           ----------------------------
    Refrigerator Perry    | 400         79         5.4
    Jerry Rice            | 180         76         4.5
    Zachary Jones         |  28         25        30.0
    Too Tall Jones        | 270         81         5.0
    Charlie Joiner        | 185         78         4.6

The data above is the "obs" matrix.  Each row in the 2D array is a data point, 
often called an "observation" or "sample".  Each column is sometimes called the 
"features" of a player. Imagine, we want to split this data set into two 
"clusters", perhaps dividing the data into linemen and receivers.  One way to 
find two "codes", one to represent each of these groups.  (I think the term 
"code" comes from communication theory, but I'm not sure.  Perhaps "classes" is 
more descriptive.)  I watched enough games (observed enough players) to make an 
educated guess as to what these codes might be:

possible code book:

                code        lb         inches     seconds
                            -----------------------------
    receiver     0         | 180          75         4.8
    lineman      1         | 260          76         5.5


So code 0 stands for a "typical" receiver and code 1 represents your typical 
lineman.  "Vector quantization" is an algorithm that calculates the distance 
between a data point and every code in the "code book" to find the closest one 
(i.e. which class is the best match). In scipy.cluster, the vq module houses 
the vector quantization tools. vq.vq(obs,code_book) returns 2 arrays -- (1) the 
index (row in the code book) of the code nearest to each data point, and (2) 
the distance that each data point is from that code. code_book is always a 2D
array.  If obs is a 2D array, each row is a separate data point.

    # note I rounded some numbers to save space
    >>> obs = array(((400, 79, 5.4),
    ...              (180, 76, 4.5),
    ...              (28,  25, 30.),
    ...              (270, 81, 5.0),
    ...              (185, 78, 4.6)))
    >>> code_book = array(((180, 75, 4.8),
    ...                    (260, 76, 5.5)))
    >>> code_id, distortion = vq.vq(obs,code_book)
    >>> code_id
    array([1, 0, 0, 1, 0])
    >>> distortion
    array([ 140.03, 1.045, 161.985, 11.192, 5.834]) 

code_id now tells what position each of the football players is most likely to 
play. Distortion is the distance, using sqrt( a^2 + b^2 + c^2), that each 
player is from the code (typical player) in their category.  Low numbers mean 
the match is good.  For example, vq tells us that Jerry Rice is receiver, and 
the low distortion means that it is pretty dang sure about that. 

                             code_id          distortion
                         ------------------------------ 
    Refrigerator Perry     1 --> lineman         140.03      
    Jerry Rice             0 --> receiver          1.04
    Zachary Jones          0 --> receiver        161.99
    Too Tall Jones         1 --> lineman          11.19
    Charlie Joiner         0 --> receiver          5.83

Most of the classifications make sense, but the distortions have some problems. 
Notably that my 1 year old son is about as likely to be a receiver as R. Perry 
is to be a lineman.  Looking at the data, it is obvious that R. Perry is a 
lineman. It isn't obvious where Zach falls because he's small (like a receiver) 
and slow (like a lineman).  So we should be quite a bit more sure about the 
fact that R. Perry is a lineman.  What's wrong?  Well, the distortion value's 
primary contribution comes from the large weight differences between these 
players and the typical players.  Even though Zach's speed is a long way from 
receiver's speed this doesn't contribute to the distortion much.  That's bad 
because the speed difference carries a lot of useful information.  For more 
accurate distortion values, we would like each feature of a player 
(weight, height, and speed) to be weighted equally.  One way of doing this is 
to "whiten" the features of both the data and the code book by dividing each 
feature by the standard deviation of that feature set.  

The standard deviation of weights and speeds across all the players are:

    #weight
    >>> stats.stdev((400,180,28,270,185))
    136.30407183939883
    #speed
    >>> stats.stdev((5.4,4.5,30.0,5.0,4.6))
    11.241885962773329

So the whitened weight and speed distortions for R. Perry from a typical 
lineman is:

    # whitened weight
    >>> abs(400-260)/136.3
    1.0271460014673512   
    # whitened speed
    >>> abs(5.4-5.5)/11.24
    0.0088967971530248789

So the whitened weight and speed distortions for Z. Jones from a typical 
receiver is:

    # whitened weight
    >>> abs(28-180)/136.3
    1.1151870873074101
    # whitened speed
    >>> (30.0-4.8)/11.24
    2.2419928825622777

It is apparent from these values that Zach's speed difference is gonna have
a large affect on the distortion now.  

The whiten() function handles the chore of whitening a data set.

    >>> wh_obs = whiten(obs)

Usually, the code book is actually calculated from a whitened set of data (via 
kmeans or some other method), and thus are automatically white.  In this case,
I specified the code_book in non-whitened units, so they'll need to be normalized
by the same standard deviation of the data set.

    # not normally needed
    >>> wh_code_book = code_book / std(obs,axis=0)

Now, rerunning vq gives:
    
    >>> code_id, distortion = vq.vq(wh_obs,wh_code_book)
    >>> code_id
    array([1, 0, 0, 1, 0])
    >>> distortion
    array([ 1.034, 0.049,  3.257 ,  0.225,  0.131])
    
                             code_id         whitened distortion
                         --------------------------------------- 
    Refrigerator Perry     1 --> lineman             1.034      
    Jerry Rice             0 --> receiver            0.049
    Zachary Jones          0 --> receiver            3.257
    Too Tall Jones         1 --> lineman             0.225
    Charlie Joiner         0 --> receiver            0.131

Now Zach's distortion is much higher than everyone elses which makes sense.

In the example above, I made an educated guess based on my knowledge of 
football of what the size and speed of a typical player in a given position 
might be.  This is information I had to supply to the clustering algorithm
before it could determine how to classify each player.  Suppose you didn't
know anything about football, and watched a football gamebut were given a list of players with there
position, weight, height, and speed.  You could use this information to 
educate yourself about the traits of a typical receiver, linebacker, lineman,
etc.  The kmeans algorithm also does this.  It takes a set of data, and
the number of positions you want to categorize 

----- Original Message ----- 
From: "Costas Malamas" <costasm@hotmail.com>
To: <scipy-user@scipy.net>
Sent: Monday, January 07, 2002 5:50 PM
Subject: [SciPy-user] Kmeans help and C source


> Hi all,
> 
> I need to use a modified K-means algorithm for a project and I was delighted 
> to discover that SciPy includes a python wrapper for a kmeans() function.
> 
> However, I am not quite following the kmeans() functionality (I am new to 
> this clustering business, so this maybe a stupid newbie question): my docs 
> tell me that kmeans should partition a dataset into k clusters.  So, I 
> expect vq.kmeans(dataset, 2) to return to me dataset split up into two 
> "equivalent" datasets.  However, anyway I feed my data into vq.kmeans() this 
> doesn't happen (e.g. I feed it a 5x4 dataset and I get back two 5x1 
> vectors).  
> My guess is that either this vq.kmeans() does something different 
> --I confess to not understanding the docstring as the observation/codebook 
> terminology has no parallel to the docs I've read-- or that I am not doing 
> something right.  Any pointers? Even some documentation on the algorithm 
> would be great help.
> 
> Secondly, as I mentioned above, I need a modified kmeans.  However, I see no 
> C/Fortran code in the src tarball or CVS that seems related to kmeans.  Is 
> the base code available?  If so, is it hackable by a SWIG newbie? (I am 
> aware of SWIG, but I have never used it for anything serious).
> 
> Any and all info will be greatly appreciated :-) --and thanks for SciPy!
> 
> 
> Costas Malamas
> 
> _________________________________________________________________
> Get your FREE download of MSN Explorer at http://explorer.msn.com/intl.asp.
> 
> _______________________________________________
> SciPy-user mailing list
> SciPy-user@scipy.net
> http://www.scipy.net/mailman/listinfo/scipy-user
> 
