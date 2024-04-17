import { Observable } from 'rxjs';

export interface BackendSocketClientInterface {
  response: Observable<any>;

  sendRequest(message: any): void;
}
