import { Observable } from 'rxjs';

export interface UserIdServiceInterface {
  userId: Observable<string | null>;

  setUserId(): Promise<void>;
}
