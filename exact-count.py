circle_memo={}
line_memo={}
line_r_memo={}
line_rr_memo={}
line_rl_memo={}
line_rrl_memo={}
line_rrll_memo={}

def exact_circle(N,dms,pts,odms):
    #counts the number of ways of putting dms dimers and pts points
    #on a circle of N sites, then choose odms dimers
    
    #Impossible configurations
    if N<0 or dms<0 or pts<0 or odms<0:
        return 0
    if 2*dms+pts > N:
        return 0
    if 2*odms > 2*dms+pts:
        return 0
    if (N,dms,pts,odms) not in circle_memo:
        #Fix some edge. On this edge, there could be a dimer present
        ##it could be chosen
        result = exact_line(N-2,dms-1,pts,odms-1)
        ##it could not be chosen
        result += exact_line_rl(N,dms-1,pts,odms)
        #There could be a dimer and a point
        ##it could be chosen
        result += 2*exact_line_r(N-2,dms-1,pts-1,odms-1)
        ##it could not be chosen
        result += 2*exact_line_rrl(N,dms-1,pts-1,odms-1)
        #There could be two points
        ##it could be chosen
        result += exact_line(N-2,dms,pts-2,odms-1)
        ##it could not be chosen
        result += exact_line_rl(N,dms,pts-2,odms)
        #There could be two dimers
        ##it could be chosen
        result += exact_line_rl(N-2,dms-2,pts,odms-1)
        ##it could not be chosen
        result += exact_line_rrll(N,dms-2,pts,odms)
        #Other configurations cannot be chosen
        #just a dimer
        result += 2*exact_line_rr(N-1,dms-1,pts,odms)
        #just a point
        result += 2*exact_line_r(N-1,dms,pts-1,odms)
        #nothing
        result += exact_line(N-2,dms,pts,odms)
        circle_memo[(N,dms,pts,odms)]= result
    return circle_memo[(N,dms,pts,odms)]

def exact_line(N,dms,pts,odms):
    #count the number of ways of choosing dms dimers and pts points
    #on a line of N sites, then choose odms dimers
    #basic configuration
    if (N,dms,pts,odms)==(0,0,0,0):
        return 1
    #impossible configurations
    if 2*dms+pts>N:
        return 0
    if 2*odms>2*dms+pts:
        return 0
    if N<0 or dms<0 or pts<0 or odms<0:
        return 0
    if (N,dms,pts,odms) not in line_memo:
        #at one end, there could be either a dimer, or a point, or nothing
        result = exact_line(N-1,dms,pts,odms)
        result += exact_line_r(N,dms,pts-1,odms)
        result += exact_line_rr(N,dms-1,pts,odms)
        line_memo[(N,dms,pts,odms)]=result
    return line_memo[(N,dms,pts,odms)]

def exact_line_r(N,dms,pts,odms):
    #dumb configurations
    if (N,dms,pts,odms)==(1,0,0,0):
        return 1
    if 2*dms+pts+1>N:
        return 0
    if 2*odms>2*dms+pts+1:
        return 0
    if N<1 or dms<0 or pts<0 or odms<0:
        return 0
    if (N,dms,pts,odms) not in line_r_memo:
        #there could be no chosen dimer at the end
        result = exact_line(N-1,dms,pts,odms)
        #there could be a chosen dimer at the end
        ##either there is a point
        result += exact_line(N-2,dms,pts-1,odms-1)
        ##or there is a dimer
        result += exact_line_r(N-2,dms-1,pts,odms-1)
        line_r_memo[(N,dms,pts,odms)]=result
    return line_r_memo[(N,dms,pts,odms)]

def exact_line_rr(N,dms,pts,odms):
    #dumb configurations
    if (N,dms,pts,odms)==(2,0,0,0):
        return 1
    if 2*dms+pts+2>N:
        return 0
    if 2*odms>2*dms+pts+2:
        return 0
    if N<2 or dms<0 or pts<0 or odms<0:
        return 0
    if (N,dms,pts,odms) not in line_rr_memo:
        #there could be no chosen dimer at the end
        result = exact_line_r(N-1,dms,pts,odms)
        #or there could be one
        result += exact_line(N-2,dms,pts,odms-1)
        line_rr_memo[(N,dms,pts,odms)]=result
    return line_rr_memo[(N,dms,pts,odms)]

def exact_line_rl(N,dms,pts,odms):
    #dumb configurations
    if (N,dms,pts,odms)==(2,0,0,0):
        return 1
    if 2*dms+pts+2>N:
        return 0
    if 2*odms>2*dms+pts+2:
        return 0
    if N<2 or dms<0 or pts<0 or odms<0:
        return 0
    #there could be no chosen dimer at one end
    result = exact_line_r(N-1,dms,pts-1,odms)
    #there could be a chosen dimer at the end
    ##either there is a point
    result += exact_line_r(N-2,dms,pts-1,odms-1)
    ##or there is a dimer
    result += exact_line_rl(N-2,dms-1,pts,odms-1)
    return result

def exact_line_rrl(N,dms,pts,odms):
    #dumb configurations
    if (N,dms,pts,odms)==(3,0,0,0):
        return 1
    if 2*dms+pts+3>N:
        return 0
    if 2*odms>2*dms+pts+3:
        return 0
    if N<3 or dms<0 or pts<0 or odms<0:
        return 0

    #there could be a chosen dimer at the end
    result = exact_line_r(N-2,dms,pts,odms-1)
    #or there could be no chosen dimer at the end
    result += exact_line_rl(N-1,dms,pts,odms)
    return result

def exact_line_rrll(N,dms,pts,odms):
    #dumb configurations
    if (N,dms,pts,odms)==(4,0,0,0):
        return 1
    if 2*dms+pts+4>N:
        return 0
    if 2*odms>2*dms+pts+4:
        return 0
    if N<4 or dms<0 or pts<0 or odms<0:
        return 0

    #there could be a chosen dimer at the end
    result = exact_line_rr(N-2,dms,pts,odms-1)
    #or there could be none
    result += exact_line_rrl(N-1,dms,pts,odms)
    return result
