"""Read Rust/C++ trusted fixture producers; compare physical matrices/actions."""
import json
import math
import pathlib
import sys
import rotations as r

def records(path):
    text=pathlib.Path(path).read_text()
    try:
        value=json.loads(text)
        return value if isinstance(value,list) else value['cases']
    except json.JSONDecodeError:
        return [json.loads(line) for line in text.splitlines() if line.strip()]

def error(a,b):
    return max(abs(x-y) for ar,br in zip(a,b) for x,y in zip(ar,br))

def verify(record):
    angles=record.get('product_angles',record['angles'])
    expected=r.from_euler(record['convention'],*angles)
    matrix=r.validate_matrix(record['matrix'])
    errors=[error(matrix,expected),error(r.from_quaternion(record['quaternion']),expected)]
    recovered=r.to_euler(matrix,record['convention'])
    errors.append(error(r.from_euler(record['convention'],*recovered.values),expected))
    axis,angle=r.to_axis_angle(matrix)
    errors.append(error(r.from_axis_angle(axis,angle),expected))
    if 'axis' in record:
        errors.append(error(r.from_axis_angle(record['axis'],record['angle']),expected))
    if 'extracted' in record:
        errors.append(error(r.from_euler(record['convention'],*record['extracted']),expected))
    if 'ck' in record:
        ck=tuple(complex(*value) for value in record['ck'])
        errors.append(error(r.ck_to_goldstein(ck),r.transpose(expected)))
    if not all(math.isfinite(e) and e<=5e-11 for e in errors):
        raise AssertionError((record['name'],errors))
    return max(errors)

results=[]
for path in sys.argv[1:]:
    data=records(path)
    if not data: raise AssertionError('empty fixture producer')
    maximum=max(verify(row) for row in data)
    results.append({'file':str(path),'cases':len(data),'max_component_error':maximum})
    # A valid proper rotation with the wrong handedness/order must still fail.
    witness=next(row for row in data if row['name'].endswith('regular') or row['name'].endswith('_0'))
    bad=dict(witness,matrix=r.transpose(witness['matrix']))
    try: verify(bad)
    except AssertionError: pass
    else: raise AssertionError('wrong transpose negative control accepted')
print(json.dumps({'status':'PASS','tolerance':5e-11,'producers':results,'comparison':'matrices and vector actions, not ambiguous quaternion signs/Euler triples','negative_control':'transposed proper rotation rejected'}))
