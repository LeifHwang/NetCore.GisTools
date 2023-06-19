//import 'bootstrap/dist/css/bootstrap.css';
import React from 'react';
import ReactDOM from 'react-dom';

import App from './App';
import registerServiceWorker from './registerServiceWorker';

const rootElement = document.getElementById('root');
ReactDOM.render(<App />, rootElement);

registerServiceWorker();

