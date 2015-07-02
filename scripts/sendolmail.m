function sendolmail(to,subject,body,attachments)
%Sends email using MS Outlook. The format of the function is
%Similar to the SENDMAIL command.

if nargin < 2
    subject = '';
end
if nargin < 3
    body = '';
end

% Create object and set parameters.
h = actxserver('outlook.Application');
mail = h.CreateItem('olMail');
mail.Subject = subject;
mail.To = to;
mail.BodyFormat = 'olFormatHTML';
mail.HTMLBody = body;

% Add attachments, if specified.
if nargin == 4
for i = 1:length(attachments)
mail.attachments.Add(attachments{i});
end
end

% Send message and release object.
mail.Send;
h.release;