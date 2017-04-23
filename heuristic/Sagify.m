intrinsic SagifyDescription(obj::.) -> MonStgElt
{Converts the description into a list that, when evaluated in Sage, corresponds
with the old list.}

case Type(obj):
    when List: return "[" cat &cat[ SagifyDescription(x) cat ",": x in obj ] cat "]";
    when SeqEnum: return "[" cat &cat[ SagifyDescription(x) cat ",": x in obj ] cat "]";
    when RngIntElt: return Sprint(obj);
    when FldRatElt: return Sprint(obj);
    when MonStgElt: return "'" cat Sprint(obj) cat "'";
end case;

end intrinsic;
