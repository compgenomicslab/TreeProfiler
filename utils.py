import re
import operator
import math
operator_dict = {
                '<':operator.lt,
                '<=':operator.le,
                '=':operator.eq,
                '!=':operator.ne,
                '>':operator.gt,
                '>=':operator.ge,
                }

def check_nan(value):
    try:
        return math.isnan (float(value))
    except ValueError:
        return False

def counter_call(node, internal_prop, leaf_prop, datatype, operator_string, right_value):
    if datatype == 'str':
        counter_props = node.props.get(internal_prop)
            
        if counter_props:
            counter_datas = counter_props.split('||')
            for counter_data in counter_datas:
                k, v = counter_data.split('--')
                if k == leaf_prop:
                    left_value = float(v)
                    return operator_dict[operator_string](left_value, float(right_value))
                else:
                    pass
        else:    
            return False
        
    else:
        return False

def call(node, prop, datatype, operator_string, right_value):
    num_operators = [ '<', '<=', '>', '>=' ] 
    if datatype == 'str':
        if operator_string in num_operators:
            return False
        elif operator_string == 'contains':
            left_value = node.props.get(prop)
            if left_value:
                return right_value in left_value 
        elif operator_string == 'in':
            left_value = right_value
            right_value = node.props.get(prop)
            if right_value:
                return left_value in right_value 
        else:
            left_value = node.props.get(prop)
            
            if left_value:
                return operator_dict[operator_string](left_value, right_value)
    
    elif datatype == 'num':
        left_value = node.props.get(prop)
        if left_value:
            return operator_dict[operator_string](float(left_value), float(right_value))
        else:
            return False

def to_code(string):
    conditional_output = []
    operators = [ '<', '<=', '>', '>=', '=', '!=', 'in', 'contains'] 
    
    r = re.compile( '|'.join( '(?:{})'.format(re.escape(o)) for o in sorted(operators, reverse=True, key=len)) )

    #codes = string.split(',')
    # code = code.replace(",", " and ") # ',' means and
    # code = code.replace(";", " or ") # ';' means or
    # code = code.replace("=", " = ")
    # code = code.replace('>', " > ")
    # code = code.replace('>=', ' >= ')

    condition_strings = string.split(',')
    for condition_string in condition_strings:
        ops = r.findall(condition_string)
        for op in ops:
            condition_string = re.sub(op, ' '+op+' ', condition_string)
            left_value, op, right_value = condition_string.split(None,2)
            conditional_output.append([left_value, op, right_value])

    return conditional_output

    