import json
with open('/Users/kevin/.gemini/antigravity/brain/63d231bc-583c-4b66-aa34-ea08742d588e/.system_generated/logs/transcript.jsonl') as f:
    for line in f:
        obj = json.loads(line)
        if obj.get('type') == 'USER_INPUT':
            print(obj.get('content'))
            break
